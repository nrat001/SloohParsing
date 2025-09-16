#!/usr/bin/env python3
"""
Summarise a directory of FITS files and (optionally) upload matched files to BHTOM.

Core features:
- Recursively scan for .fit/.fits/.fts under a root folder (e.g., data/NotUploaded).
- Parse RA/Dec from filename prefix: HHMMSS[p|m]DDMMSS_...
- Read FITS header keys: DATE-OBS, JD, FILTER, EXPTIME, OBSERVER.
- Match filename coords to lookup CSV (name,ra_hours_decimal,dec_deg_decimal).
- Write CSV/XLSX for files that REMAIN in the scanned folder at end of processing
  (i.e., files that were not moved due to a successful upload).
- (Optional) Upload matched files to BHTOM; on success (and not --dry-run),
  move file from NotUploaded/ -> Uploaded/ (sibling folder).

OBSERVATORY/FILTER logic:
- Observatory name for upload is determined by:
  1) CLI override --observatory if provided; else
  2) OBSERVER header -> mapped name via OBSERVER_TO_OBSNAME (unknowns are skipped).
- Filter for upload is determined by:
  1) CLI override --filter if provided; else
  2) OBSERVER header -> mapped filter via map_filter_from_observer:
     - Australia One, Chile Two, Canary Two -> 'GaiaSP/UBVRI'
     - Canary One -> 'GaiaSP/any'
     Unknowns are skipped.

CLI:
  python3 sloohParser.py FITS_DIR LOOKUP_CSV
      [--csv OUT.csv] [--xlsx OUT.xlsx]
      [--tol-arcmin N]
      [--upload] [--dry-run]
      [--filter FILTER_NAME] [--observatory OBSNAME]
      [--token-env ENVVAR]

Typical:
  export BHTOM_TOKEN='<your token>'
  python3 sloohParser.py ./data/NotUploaded ./lookup.csv \
    --xlsx ./fits_summary.xlsx \
    --upload --dry-run
"""

import argparse
import csv
import math
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests
from astropy.io import fits  # headers only
import pandas as pd

# ------------------------- Config -------------------------

# Filename pattern for RA/Dec prefix: HHMMSS [p|m] DDMMSS _
FILENAME_REGEX = re.compile(
    r'(?P<rah>\d{2})(?P<ram>\d{2})(?P<ras>\d{2})(?P<sign>[pmPM])(?P<dd>\d{2})(?P<dm>\d{2})(?P<ds>\d{2})_'
)

MATCH_TOL_ARCMIN_DEFAULT = 10.0

# Output columns (Target name first). We will write rows for files that remain in place after processing.
CSV_COLUMNS = [
    "Target name",
    "Filename",
    "DATE-OBS",
    "JD",
    "Filter",
    "EXPTIME",
    "OBSERVER",
    "RA_hours",
    "RA_min",
    "RA_sec",
    "DEC_deg",
    "DEC_min",
    "DEC_sec",
    "RA_hours_decimal",
    "DEC_deg_decimal",
    "Match_sep_arcmin",
    "Cat_RA_hours_decimal",
    "Cat_DEC_deg_decimal",
]

# BHTOM upload service endpoint
UPLOAD_URL = "https://uploadsvc2.astrolabs.pl/upload/"

# Known OBSERVER header -> BHTOM observatory_name
OBSERVER_TO_OBSNAME: Dict[str, str] = {
    "Slooh.com Chile Two Wide-Field Telescope": "Slooh-CL_FLI-KAF-16803",
    "Slooh.com Australia One Half Meter Telescope": "Slooh-AU_FLI-KAF-16803",
    # Add Canary mappings as you learn them, e.g.:
    # "Slooh.com Canary One ...": "Slooh-??_....",
    # "Slooh.com Canary Two ...": "Slooh-??_....",
}

# ---------------------- Helpers (coords) ----------------------

def parse_filename_coords(fname: str) -> Optional[Dict[str, int]]:
    m = FILENAME_REGEX.search(fname)
    if not m:
        return None
    rah = int(m.group("rah"))
    ram = int(m.group("ram"))
    ras = int(m.group("ras"))
    sign = -1 if m.group("sign").lower() == "m" else 1
    dd = int(m.group("dd"))
    dm = int(m.group("dm"))
    ds = int(m.group("ds"))
    return {"rah": rah, "ram": ram, "ras": ras, "sign": sign, "dd": dd, "dm": dm, "ds": ds}

def hms_to_hours(h: int, m: int, s: int) -> float:
    return h + m/60.0 + s/3600.0

def dms_to_degrees(sign: int, d: int, m: int, s: int) -> float:
    return float(sign) * (d + m/60.0 + s/3600.0)

def ang_sep_arcmin(ra1_hours: float, dec1_deg: float, ra2_hours: float, dec2_deg: float) -> float:
    """Great-circle separation in arcminutes between (ra1,dec1) and (ra2,dec2)."""
    ra1 = math.radians(ra1_hours * 15.0)
    ra2 = math.radians(ra2_hours * 15.0)
    dec1 = math.radians(dec1_deg)
    dec2 = math.radians(dec2_deg)
    d_ra = ra2 - ra1
    d_dec = dec2 - dec1
    sin2 = math.sin(d_dec/2.0)**2 + math.cos(dec1)*math.cos(dec2)*math.sin(d_ra/2.0)**2
    if sin2 < 0.0:
        sin2 = 0.0
    elif sin2 > 1.0:
        sin2 = 1.0
    d_sigma = 2.0 * math.asin(math.sqrt(sin2))  # radians
    return math.degrees(d_sigma) * 60.0

# ---------------------- FITS header helpers ----------------------

def read_header_cards(path: Path) -> Dict[str, str]:
    """Return selected header keys (missing -> '')."""
    keys = ["DATE-OBS", "JD", "FILTER", "EXPTIME", "OBSERVER"]
    out: Dict[str, str] = {k: "" for k in keys}
    try:
        hdr = fits.getheader(path, 0)
    except Exception:
        return out

    def get_any(*cands: str) -> str:
        for k in cands:
            if k in hdr and hdr[k] is not None:
                return str(hdr[k]).strip()
        return ""

    out["DATE-OBS"] = get_any("DATE-OBS", "DATEOBS", "DATE")
    out["JD"]       = get_any("JD", "MJD-OBS", "MJD")
    out["FILTER"]   = get_any("FILTER", "FILTER1", "FILTER2")
    out["EXPTIME"]  = get_any("EXPTIME", "EXPOSURE", "ITIME")
    out["OBSERVER"] = get_any("OBSERVER", "OBSERVER1", "AUTHOR")
    return out

# ---------------------- Catalog & matching ----------------------

def load_lookup(lookup_csv: Path) -> List[Tuple[str, float, float]]:
    """Load catalog rows as (name, ra_hours, dec_deg)."""
    rows: List[Tuple[str, float, float]] = []
    with lookup_csv.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            name = str(row["name"]).strip()
            ra_h = float(row["ra_hours_decimal"])
            dec_d = float(row["dec_deg_decimal"])
            rows.append((name, ra_h, dec_d))
    return rows

def match_target(
    ra_hours: float,
    dec_deg: float,
    catalog: List[Tuple[str, float, float]],
    tol_arcmin: float,
) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[float]]:
    """Return (name, sep_arcmin, cat_ra_hours, cat_dec_deg) if best match is within tol; else all None."""
    best_name: Optional[str] = None
    best_sep: Optional[float] = None
    best_ra: Optional[float] = None
    best_dec: Optional[float] = None

    for name, cat_ra_h, cat_dec_d in catalog:
        sep = ang_sep_arcmin(ra_hours, dec_deg, cat_ra_h, cat_dec_d)
        if (best_sep is None) or (sep < best_sep):
            best_sep = sep
            best_name = name
            best_ra = cat_ra_h
            best_dec = cat_dec_d

    if (best_sep is not None) and (best_sep <= tol_arcmin):
        return best_name, best_sep, best_ra, best_dec
    return None, None, None, None

# ---------------------- Uploading ----------------------

def map_observatory_name(observer_header: str) -> Optional[str]:
    """Return BHTOM observatory_name for a given OBSERVER header, or None if unknown."""
    if not observer_header:
        return None
    key = observer_header.strip()
    return OBSERVER_TO_OBSNAME.get(key)

def map_filter_from_observer(observer_header: str) -> Optional[str]:
    """
    Returns the recommended filter for uploads based on the OBSERVER header.
    Rules:
      - Australia One, Chile Two, Canary Two -> 'GaiaSP/UBVRI'
      - Canary One -> 'GaiaSP/any'
    """
    if not observer_header:
        return None
    s = observer_header.strip().lower()

    # Specific known strings
    if s == "slooh.com australia one half meter telescope":
        return "GaiaSP/UBVRI"
    if s == "slooh.com chile two wide-field telescope":
        return "GaiaSP/UBVRI"

    # Canary telescopes via keywords
    if "canary" in s:
        if "one" in s:
            return "GaiaSP/any"
        if "two" in s:
            return "GaiaSP/UBVRI"

    # Unknown
    return None

def upload_fits(
    file_path: Path,
    target_name: str,
    observatory_name: str,
    filter_name: str,
    token: str,
    dry_run: bool = False,
) -> dict:
    """
    Upload a single FITS file to BHTOM upload service.
    Returns parsed JSON (dict). Raises on non-2xx status.
    """
    headers = {"Authorization": f"Token {token}"}
    data = {
        "target": target_name,
        "observatory": observatory_name,
        "filter": filter_name,
        "data_product_type": "fits_file",
    }
    if dry_run:
        data["dry_run"] = "True"  # API expects string 'True'

    with open(file_path, "rb") as fh:
        files = {"files": (file_path.name, fh, "application/fits")}
        resp = requests.post(UPLOAD_URL, headers=headers, data=data, files=files, timeout=180)
    resp.raise_for_status()
    return resp.json()

# ---------------------- File discovery & writers ----------------------

def find_fits_files(root: Path) -> List[Path]:
    """Find .fit/.fits/.fts files recursively under root."""
    exts = (".fit", ".fits", ".fts")
    return [p for p in root.rglob("*") if p.suffix.lower() in exts and p.is_file()]

def write_csv(rows: List[Dict[str, object]], out_path: Path) -> None:
    with out_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

def write_xlsx(rows: List[Dict[str, object]], out_path: Path) -> None:
    df = pd.DataFrame(rows, columns=CSV_COLUMNS)
    engine = None
    try:
        import openpyxl  # noqa: F401
        engine = "openpyxl"
    except Exception:
        try:
            import xlsxwriter  # noqa: F401
            engine = "xlsxwriter"
        except Exception:
            engine = None
    if engine is None:
        raise RuntimeError("No Excel writer engine found. Install either 'openpyxl' or 'xlsxwriter'.")
    with pd.ExcelWriter(out_path, engine=engine) as writer:
        df.to_excel(writer, sheet_name="summary", index=False)

# --------------------------- Main ---------------------------

def process_file(
    path: Path,
    catalog: List[Tuple[str, float, float]],
    tol_arcmin: float,
) -> Optional[Dict[str, object]]:
    """Build a result row for a single file (no uploading/moving here)."""
    parts = parse_filename_coords(path.name)
    if parts is None:
        return None

    # Decimal coords from filename
    ra_hours_dec = hms_to_hours(parts["rah"], parts["ram"], parts["ras"])
    dec_deg_dec  = dms_to_degrees(parts["sign"], parts["dd"], parts["dm"], parts["ds"])

    target_name, best_sep, cat_ra, cat_dec = match_target(ra_hours_dec, dec_deg_dec, catalog, tol_arcmin)
    hdr = read_header_cards(path)

    out_target = target_name if target_name is not None else "NO MATCH"

    row = {
        "Target name": out_target,
        "Filename": path.name,
        "DATE-OBS": hdr["DATE-OBS"],
        "JD": hdr["JD"],
        "Filter": hdr["FILTER"],
        "EXPTIME": hdr["EXPTIME"],
        "OBSERVER": hdr["OBSERVER"],
        "RA_hours": parts["rah"],
        "RA_min": parts["ram"],
        "RA_sec": parts["ras"],
        "DEC_deg": parts["sign"] * parts["dd"],  # signed integer
        "DEC_min": parts["dm"],
        "DEC_sec": parts["ds"],
        "RA_hours_decimal": round(ra_hours_dec, 6),
        "DEC_deg_decimal": round(dec_deg_dec, 6),
        "Match_sep_arcmin": round(best_sep, 3) if best_sep is not None else "",
        "Cat_RA_hours_decimal": round(cat_ra, 6) if cat_ra is not None else "",
        "Cat_DEC_deg_decimal": round(cat_dec, 6) if cat_dec is not None else "",
    }
    return row

def main() -> None:
    ap = argparse.ArgumentParser(description="Summarise FITS files and (optionally) upload matched files to BHTOM.")
    ap.add_argument("fits_dir", type=Path, help="Root directory containing FITS files (e.g., data/NotUploaded)")
    ap.add_argument("lookup_csv", type=Path, help="CSV with columns: name,ra_hours_decimal,dec_deg_decimal")

    ap.add_argument("--csv", type=Path, default=Path("fits_summary.csv"),
                    help="Output CSV path (default: fits_summary.csv)")
    ap.add_argument("--xlsx", type=Path, default=Path("fits_summary.xlsx"),
                    help="Output Excel path (default: fits_summary.xlsx)")
    ap.add_argument("--tol-arcmin", type=float, default=MATCH_TOL_ARCMIN_DEFAULT,
                    help=f"Match tolerance in arcmin (default: {MATCH_TOL_ARCMIN_DEFAULT})")

    # Upload options
    ap.add_argument("--upload", action="store_true",
                    help="If set, upload each matched FITS to BHTOM.")
    ap.add_argument("--dry-run", action="store_true",
                    help="If set with --upload, perform a dry-run upload (server validates but does not persist).")
    ap.add_argument("--observatory", type=str, default="",
                    help="Observatory name for BHTOM (optional; overrides OBSERVER->observatory mapping).")
    ap.add_argument("--filter", type=str, default="",
                    help="Filter name recognized by BHTOM (optional; if omitted, derived from OBSERVER when possible).")
    ap.add_argument("--token-env", type=str, default="BHTOM_TOKEN",
                    help="Environment variable holding the BHTOM API token (default: BHTOM_TOKEN).")

    args = ap.parse_args()

    # Load catalog
    catalog = load_lookup(args.lookup_csv)

    # Find files to process (these are candidates to remain if not moved)
    files = sorted(find_fits_files(args.fits_dir))

    # Summary counters
    total_scanned = len(files)
    parsed_ok = 0
    unmatched_count = 0
    upload_attempted = 0
    upload_success = 0
    upload_dryrun_success = 0
    upload_failed = 0
    skipped_unknown_observer = 0
    skipped_no_match = 0
    skipped_no_filter = 0

    unknown_observer_values: set[str] = set()

    # Rows to write at the end (only files that REMAIN in the scanned dir)
    rows_remaining: List[Dict[str, object]] = []

    # Prepare upload destination for moves (if we detect a NotUploaded parent)
    # If file path looks like .../NotUploaded/<file>, move to sibling .../Uploaded/
    def destination_for(f: Path) -> Path:
        parent = f.parent
        if parent.name == "NotUploaded":
            uploaded_dir = parent.parent / "Uploaded"
        else:
            # Default: create an "Uploaded" sibling of the current folder
            uploaded_dir = parent / "Uploaded"
        uploaded_dir.mkdir(parents=True, exist_ok=True)
        return uploaded_dir / f.name

    # Token (if uploading)
    token = ""
    if args.upload:
        token = os.environ.get(args.token_env, "").strip()
        if not token:
            raise RuntimeError(f"Env var {args.token_env} is empty or not set; cannot upload.")
        # --filter is optional; if not provided we map per-file from OBSERVER.

    for f in files:
        try:
            row = process_file(f, catalog, float(args.tol_arcmin))
            if row is None:
                # Could not parse filename; keep it in the "remaining" outputs for visibility
                rows_remaining.append({
                    "Target name": "NO MATCH",
                    "Filename": f.name,
                    "DATE-OBS": "",
                    "JD": "",
                    "Filter": "",
                    "EXPTIME": "",
                    "OBSERVER": "",
                    "RA_hours": "",
                    "RA_min": "",
                    "RA_sec": "",
                    "DEC_deg": "",
                    "DEC_min": "",
                    "DEC_sec": "",
                    "RA_hours_decimal": "",
                    "DEC_deg_decimal": "",
                    "Match_sep_arcmin": "",
                    "Cat_RA_hours_decimal": "",
                    "Cat_DEC_deg_decimal": "",
                })
                continue

            parsed_ok += 1
            target_name = str(row.get("Target name", "NO MATCH"))

            if args.upload:
                if target_name == "NO MATCH":
                    skipped_no_match += 1
                    # Keep it in remaining outputs
                    rows_remaining.append(row)
                    print(f"[SKIP] {f.name}: no target match; not uploading.")
                else:
                    # Decide observatory_name: manual override or map from OBSERVER
                    observatory_name = args.observatory.strip() if args.observatory.strip() \
                        else map_observatory_name(str(row.get("OBSERVER", "")))

                    if not observatory_name:
                        skipped_unknown_observer += 1
                        unknown_observer_values.add(str(row.get("OBSERVER", "")) or "<empty>")
                        rows_remaining.append(row)
                        print(f"[SKIP] {f.name}: unknown OBSERVER='{row.get('OBSERVER','')}'. "
                              f"Add mapping in OBSERVER_TO_OBSNAME or pass --observatory to override.")
                    else:
                        # Decide filter: CLI override wins; otherwise map from OBSERVER
                        upload_filter = args.filter.strip() if args.filter.strip() \
                                        else map_filter_from_observer(str(row.get("OBSERVER", "")))

                        if not upload_filter:
                            skipped_no_filter += 1
                            rows_remaining.append(row)
                            print(f"[SKIP] {f.name}: no filter mapping for OBSERVER='{row.get('OBSERVER','')}'. "
                                  f"Pass --filter to override.")
                        else:
                            upload_attempted += 1
                            try:
                                result = upload_fits(
                                    file_path=f,
                                    target_name=target_name,
                                    observatory_name=observatory_name,
                                    filter_name=upload_filter,
                                    token=token,
                                    dry_run=bool(args.dry_run),
                                )
                                success = True
                                if isinstance(result, dict) and "Success" in result:
                                    success = True
                                if args.dry_run:
                                    upload_dryrun_success += 1
                                    # Dry-run: do NOT move file; keep it in remaining outputs
                                    rows_remaining.append(row)
                                    print(f"[OK][DRY-RUN] {f.name} ({observatory_name}, {upload_filter}): {result}")
                                else:
                                    if success:
                                        upload_success += 1
                                        # Move the file to Uploaded/
                                        dest = destination_for(f)
                                        try:
                                            f.rename(dest)
                                            print(f"[OK][MOVED] {f.name} -> {dest} (obs={observatory_name}, filter={upload_filter})")
                                        except Exception as move_err:
                                            print(f"[WARN] Uploaded but failed to move {f.name}: {move_err}")
                                            # If move failed, the file still remains here; include it in outputs
                                            rows_remaining.append(row)
                                    else:
                                        upload_failed += 1
                                        rows_remaining.append(row)
                                        print(f"[ERR] Upload response indicates failure for {f.name}: {result}")
                            except requests.HTTPError as e:
                                upload_failed += 1
                                rows_remaining.append(row)
                                print(f"[ERR] Upload failed for {f.name}: HTTP {e.response.status_code} {e.response.text}")
                            except Exception as e:
                                upload_failed += 1
                                rows_remaining.append(row)
                                print(f"[ERR] Upload failed for {f.name}: {e}")
            else:
                # Not uploading: keep the original behavior (list all found)
                rows_remaining.append(row)

            if target_name == "NO MATCH":
                unmatched_count += 1

        except Exception as e:
            # Any unexpected error: keep file in outputs for visibility
            rows_remaining.append({
                "Target name": "NO MATCH",
                "Filename": f.name,
                "DATE-OBS": "",
                "JD": "",
                "Filter": "",
                "EXPTIME": "",
                "OBSERVER": "",
                "RA_hours": "",
                "RA_min": "",
                "RA_sec": "",
                "DEC_deg": "",
                "DEC_min": "",
                "DEC_sec": "",
                "RA_hours_decimal": "",
                "DEC_deg_decimal": "",
                "Match_sep_arcmin": "",
                "Cat_RA_hours_decimal": "",
                "Cat_DEC_deg_decimal": "",
            })
            print(f"[ERR] Processing failed for {f.name}: {e}")

    # Write remaining files summary (what still resides in the scanned directory)
    write_csv(rows_remaining, args.csv)
    try:
        write_xlsx(rows_remaining, args.xlsx)
    except Exception as e:
        print(f"WARNING: XLSX not written: {e}")

    # Final summary
    moved_count = upload_success  # each successful non-dry-run upload triggers a move
    remaining_count = len(rows_remaining)

    print("\n=== Summary ===")
    print(f"Scanned files:                {total_scanned}")
    print(f"Parsed OK:                    {parsed_ok}")
    print(f"Unmatched (NO MATCH):         {unmatched_count}")
    if args.upload:
        print(f"Upload attempted:             {upload_attempted}")
        if args.dry_run:
            print(f"  Dry-run successes:          {upload_dryrun_success}")
            print(f"  Files moved:                0 (dry-run)")
        else:
            print(f"  Upload successes:           {upload_success}")
            print(f"  Upload failed:              {upload_failed}")
            print(f"  Files moved to 'Uploaded':  {moved_count}")
        print(f"Skipped (unknown OBSERVER):   {skipped_unknown_observer}")
        print(f"Skipped (no filter mapping):  {skipped_no_filter}")
        if unknown_observer_values:
            print("Unknown OBSERVER values seen:")
            for val in sorted(unknown_observer_values):
                print(f"  - {val}")
    print(f"Remaining in scanned dir:     {remaining_count}")
    print(f"Wrote CSV:                    {args.csv}")
    print(f"Wrote XLSX:                   {args.xlsx}")

if __name__ == "__main__":
    main()
