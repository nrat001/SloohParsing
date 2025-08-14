<<<<<<< HEAD
# FITS Summary Parser

A small command-line tool that walks a directory of astronomical FITS images, extracts key metadata from **FITS headers** and **filenames**, optionally **matches** each image to a target from a lookup table, and writes a tidy **CSV** and/or **Excel (XLSX)** summary.

## What it does

* Recursively finds `.fit`, `.fits`, `.fts` files in a given directory.
* Reads selected header cards (without loading image data): `DATE-OBS`, `JD`, `FILTER`, `EXPTIME`, `OBSERVER`.
* Parses **RA/Dec** from the filename prefix formatted as:

  ```
  HHMMSS[p|m]DDMMSS_...
  ```

  where `p` = north (`+`), `m` = south (`-`) for the declination sign.
* Writes both the **H/M/S / D/M/S** components and the **decimal** coordinates (RA in hours, Dec in degrees).
* Cross-matches each file’s coordinates to a **lookup CSV** of targets; on success outputs the **Target name**, the **catalog coordinates**, and the **angular separation** (arcminutes).
* If no match is within tolerance, sets **Target name** to **`NO MATCH`**.
* Outputs results in a consistent column order to **CSV** and/or **XLSX**.

## Quick start

```bash
# python 3.8+ recommended
pip install astropy pandas openpyxl   # openpyxl for XLSX output

python sloohParser.py /path/to/fits_dir /path/to/lookup.csv \
  --xlsx /path/to/fits_summary.xlsx \
  --csv  /path/to/fits_summary.csv \
  --tol-arcmin 10
```

* `--xlsx` and/or `--csv` are optional; if omitted the script defaults to:

  * `fits_summary.xlsx` and/or `fits_summary.csv` in the current directory (depending on which flags you provide).
* `--tol-arcmin` sets the maximum separation for an acceptable match (default **10 arcmin**).

## Lookup CSV format

The lookup table must have a header row with **exactly** these columns:

```csv
name,ra_hours_decimal,dec_deg_decimal
```

* **`name`** — human-readable target name.
* **`ra_hours_decimal`** — RA in **decimal hours** (e.g. `23.4233889`).
* **`dec_deg_decimal`** — Dec in **decimal degrees** (e.g. `-38.4469722`).

Example:

```csv
name,ra_hours_decimal,dec_deg_decimal
Gaia21fji,21.0661389,50.2445
IRAS 23226-3843,23.4233889,-38.4469722
1RXS_J173546.9-302859,17.5964167,-30.4828889
Gaia24csq,19.1140278,-26.9395833
```

## Filename convention

The parser expects filenames whose **first token** encodes coordinates:

```
<RA:HHMMSS><sign:p|m><DEC:DDMMSS>_<YYYYMMDD>_...
```

Example: `232524m382649_20250806_161626_0_yfcyd0_e_cal.fit`

* RA → `23h 25m 24s`
* Dec → `-38° 26' 49"` (because `m` indicates negative declination)
* The underscore `_` after `DDMMSS` is required by the default regex.

> Have variants? You can relax the regex in `FILENAME_REGEX`.

## Output columns (order)

1. **Target name** (`"NO MATCH"` if none within tolerance)
2. Filename
3. DATE-OBS
4. JD
5. Filter
6. EXPTIME
7. OBSERVER
8. RA\_hours (int, from filename)
9. RA\_min (int)
10. RA\_sec (int)
11. DEC\_deg (signed int, from filename sign)
12. DEC\_min (int)
13. DEC\_sec (int)
14. RA\_hours\_decimal (float, from filename)
15. DEC\_deg\_decimal (float, from filename)
16. Match\_sep\_arcmin (float; nearest separation)
17. Cat\_RA\_hours\_decimal (float; from lookup)
18. Cat\_DEC\_deg\_decimal (float; from lookup)

## Matching details

* The script computes great-circle separations with the haversine formula (no `astropy.units` needed).
* A target is considered a match if the nearest catalog entry is within `--tol-arcmin` arcminutes.
* When unmatched:

  * **Target name** is `NO MATCH`.
  * **Match\_sep\_arcmin** is left blank by default (easy to change if you want diagnostics for near-misses).

## Requirements

* Python **3.8+**
* Packages:

  * `astropy` (for reading FITS headers only)
  * `pandas` (for DataFrame/XLSX writing)
  * **One** Excel engine if you need XLSX:

    * `openpyxl` (recommended) or `xlsxwriter`

Install:

```bash
pip install astropy pandas openpyxl
# or, if you prefer:
# pip install xlsxwriter
```

## Examples

Write **only** an Excel file:

```bash
python sloohParser.py ./data ./lookup.csv --xlsx ./fits_summary.xlsx
```

Write **only** a CSV:

```bash
python sloohParser.py ./data ./lookup.csv --csv ./fits_summary.csv
```

Write **both**:

```bash
python sloohParser.py ./data ./lookup.csv \
  --csv ./fits_summary.csv \
  --xlsx ./fits_summary.xlsx
```

Increase the match tolerance to 20 arcmin:

```bash
python sloohParser.py ./data ./lookup.csv --xlsx ./fits_summary.xlsx --tol-arcmin 20
```

## Troubleshooting

* **Lookup KeyError (e.g., `'name'` or `'ra_hours_decimal'`)**
  Ensure your lookup header row matches exactly:
  `name,ra_hours_decimal,dec_deg_decimal`.

* **“NO MATCH” for many rows**

  * Increase `--tol-arcmin` (e.g., 20–30) to diagnose.
  * Verify your lookup RA is in **hours**, not degrees.
  * Check your filenames match the expected coordinate pattern.

* **NumPy `longdouble` warning**
  Benign on some platforms. You can suppress with:
  `PYTHONWARNINGS="ignore"`.

* **XLSX engine not found**
  Install `openpyxl` (`pip install openpyxl`) or use `--csv` only.

## License

MIT (recommendation). Update to your preferred license.

---

If you want the tool to **always print the nearest separation** even on “NO MATCH”, or to **sort** the output by separation for QA, those are small tweaks—open an issue or PR!
=======
nrat001@SC448046:~/SloohParsing/src$ python3 sloohParser.py /home/nrat001/SloohParsing/data     /home/nrat001/SloohParsing/config/targetlist.csv     --csv /home/nrat001/SloohParsing/output/fits_summary.csv 
>>>>>>> 673d5fdfce56caefa427e609ecd5860af8f38a31
