#!/usr/bin/env python3
"""
Build simplified ICESat-2 track geometries as LineStrings.

Usage:
    python build_is2_tracks_geometries.py /path/to/ATL03_files

- Reads each .h5 ATL03/ATL07 file
- Extracts lat/lon from a strong beam
- Downsamples to ~10 points (configurable)
- Outputs a CSV with WKT geometry in the same folder
"""

import os
import sys
import h5py
import numpy as np
import pandas as pd
from shapely.geometry import LineString
from datetime import datetime
from tqdm import tqdm

# Optional helper
sys.path.append(".")
try:
    from src.utils.Is2_utils import get_beam_strengths
except Exception:
    get_beam_strengths = None


def pick_strong_beam(h5path):
    """Try to pick a strong beam safely (uses helper if available)."""
    if get_beam_strengths:
        try:
            strong_beams, _ = get_beam_strengths(h5path)
            if strong_beams:
                return strong_beams[0]
        except Exception as e:
            print(f"⚠️ get_beam_strengths failed for {os.path.basename(h5path)}: {e}")

    # Fallback: inspect HDF5 groups
    try:
        with h5py.File(h5path, "r") as f:
            beams = [k for k in f.keys() if k.startswith("gt")]
    except Exception as e:
        print(f"⚠️ Could not open {os.path.basename(h5path)}: {e}")
        return None

    for b in beams:
        if b.endswith("l"):
            return b
    return beams[0] if beams else None


def extract_track(h5path):
    """Extract a simplified LineString from a strong beam."""
    fname = os.path.basename(h5path)
    parts = fname.replace(".h5", "").split("_")
    product = next((p for p in parts if p.startswith("ATL")), None)
    date_str = parts[1] if len(parts) > 1 else None
    date = None
    try:
        date = datetime.strptime(date_str, "%Y%m%d%H%M%S")
    except Exception:
        pass

    beam = pick_strong_beam(h5path)
    if not beam:
        return None

    try:
        with h5py.File(h5path, "r") as f:
            if "ATL03" in product:
                lat = f[f"{beam}/heights/lat_ph"][::200]
                lon = f[f"{beam}/heights/lon_ph"][::200]
            elif "ATL07" in product:
                lat = f[f"{beam}/sea_ice_segments/latitude"][::50]
                lon = f[f"{beam}/sea_ice_segments/longitude"][::50]
            else:
                return None
    except Exception as e:
        print(f"⚠️ Could not read {fname}: {e}")
        return None

    # Clean and downsample
    m = np.isfinite(lat) & np.isfinite(lon)
    lat, lon = lat[m], lon[m]
    if len(lat) < 2:
        return None

    lon = ((lon + 180) % 360) - 180  # normalize

    line = LineString(zip(lon, lat))
    return {
        "filename": fname,
        "product": product,
        "date": date,
        "geometry_wkt": line.wkt,
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: python build_is2_tracks_geometries.py /path/to/ATL03_files")
        sys.exit(1)

    folder = sys.argv[1]
    if not os.path.isdir(folder):
        print(f"Error: {folder} is not a valid directory.")
        sys.exit(1)

    files = [os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".h5")]
    records = []

    print(f"Processing {len(files)} ICESat-2 files in {folder}...")
    for fpath in tqdm(files):
        rec = extract_track(fpath, max_points=10)
        if rec:
            records.append(rec)

    if not records:
        print("No valid tracks found.")
        sys.exit(0)

    df = pd.DataFrame(records)
    output_path = os.path.join(folder, "ICESat2_tracks.csv")
    df.to_csv(output_path, index=False)
    print(f"✅ Saved {len(df)} tracks to {output_path}")


if __name__ == "__main__":
    main()
