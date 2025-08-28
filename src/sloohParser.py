#!/usr/bin/env python3
"""
Summarise a directory of FITS files into CSV and XLSX.

Outputs per file:
- Target name (first column) — "NO MATCH" if no catalog match within tolerance
- Filename, header fields (DATE-OBS, JD, Filter, EXPTIME, OBSERVER)
- RA/Dec components parsed from filename (ints)
- Decimal RA (hours) and Dec (deg) from filename
- Match separation (arcmin) and matched catalog coords

Filename prefix pattern parsed:
  <RA:HHMMSS><sign:p|m><DEC:DDMMSS>_<YYYYMMDD>_...
Example:
  232524m382649_20250806_161626_0_yfcyd0_e_cal.fit

Lookup CSV columns (header row required):
  name,ra_hours_decimal,dec_deg_decimal
e.g.
  IRAS 23226-3843,23.4233889,-38.4469722
"""

import argparse
import csv
import math
import re
from pathlib import Path

from astropy.io import fits  # for headers only
import pandas as pd          # for XLSX output

# ---------- Config ----------

FILENAME_REGEX = re.compile(
    r'(?P<rah>\d{2})(?P<ram>\d{2})(?P<ras>\d{2})(?P<sign>[pmPM])(?P<dd>\d{2})(?P<dm>\d{2})(?P<ds>\d{2})_'
)

MATCH_TOL_ARCMIN_DEFAULT = 1.0

# Target name first in output order
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

# ---------- Helpers ----------

def parse_filename_coords(fname):
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

def hms_to_hours(h, m, s):
    return h + m/60.0 + s/3600.0

def dms_to_degrees(sign, d, m, s):
    return float(sign) * (d + m/60.0 + s/3600.0)

def ang_sep_arcmin(ra1_hours, dec1_deg, ra2_hours, dec2_deg):
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

def read_header_cards(path):
    """Return selected header keys (missing -> '')."""
    keys = ["DATE-OBS", "JD", "FILTER", "EXPTIME", "OBSERVER"]
    out = {k: "" for k in keys}
    try:
        hdr = fits.getheader(path, 0)
    except Exception:
        return out

    def get_any(*cands):
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

def load_lookup(lookup_csv):
    """Load catalog rows as (name, ra_hours, dec_deg)."""
    rows = []
    with lookup_csv.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            name = str(row["name"]).strip()
            ra_h = float(row["ra_hours_decimal"])
            dec_d = float(row["dec_deg_decimal"])
            rows.append((name, ra_h, dec_d))
    return rows

def match_target(ra_hours, dec_deg, catalog, tol_arcmin):
    """Return (name, sep_arcmin, cat_ra_hours, cat_dec_deg) if best match is within tol; else all None."""
    best_name = None
    best_sep = None
    best_ra = None
    best_dec = None

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

def process_file(path, catalog, tol_arcmin):
    parts = parse_filename_coords(path.name)
    if parts is None:
        return None

    # Decimal coords from filename
    ra_hours_dec = hms_to_hours(parts["rah"], parts["ram"], parts["ras"])
    dec_deg_dec  = dms_to_degrees(parts["sign"], parts["dd"], parts["dm"], parts["ds"])

    target_name, best_sep, cat_ra, cat_dec = match_target(ra_hours_dec, dec_deg_dec, catalog, tol_arcmin)
    hdr = read_header_cards(path)

    # If no match within tolerance, show "NO MATCH"
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

def find_fits_files(root):
    exts = (".fit", ".fits", ".fts")
    return [p for p in root.rglob("*") if p.suffix.lower() in exts and p.is_file()]

def write_csv(rows, out_path):
    with out_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=CSV_COLUMNS)
    # Ensure we always write columns in the specified order
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

def write_xlsx(rows, out_path):
    df = pd.DataFrame(rows, columns=CSV_COLUMNS)
    # Pick an available engine
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

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description="Summarise FITS files to CSV and XLSX.")
    ap.add_argument("fits_dir", type=Path, help="Root directory containing FITS files")
    ap.add_argument("lookup_csv", type=Path, help="CSV: name,ra_hours_decimal,dec_deg_decimal")
    ap.add_argument("--csv", type=Path, default=Path("fits_summary.csv"),
                    help="Output CSV path (default: fits_summary.csv)")
    ap.add_argument("--xlsx", type=Path, default=Path("fits_summary.xlsx"),
                    help="Output Excel path (default: fits_summary.xlsx)")
    ap.add_argument("--tol-arcmin", type=float, default=MATCH_TOL_ARCMIN_DEFAULT,
                    help="Match tolerance in arcmin (default: %.1f)" % MATCH_TOL_ARCMIN_DEFAULT)
    args = ap.parse_args()

    catalog = load_lookup(args.lookup_csv)
    files = find_fits_files(args.fits_dir)

    rows = []
    n_bad = 0
    n_unmatched = 0

    for f in sorted(files):
        try:
            row = process_file(f, catalog, float(args.tol_arcmin))
            if row is None:
                n_bad += 1
                continue
            if row["Target name"] == "NO MATCH":
                n_unmatched += 1
            rows.append(row)
        except Exception:
            n_bad += 1
            continue

    write_csv(rows, args.csv)
    try:
        write_xlsx(rows, args.xlsx)
    except Exception as e:
        print(f"WARNING: XLSX not written: {e}")

    print(f"Wrote CSV:  {args.csv}  (rows: {len(rows)})")
    print(f"Wrote XLSX: {args.xlsx}  (rows: {len(rows)})")
    if n_bad:
        print(f"Skipped {n_bad} file(s) due to parse/read issues.")
    if n_unmatched:
        print(f"{n_unmatched} file(s) had NO MATCH within tolerance.")

if __name__ == "__main__":
    main()
