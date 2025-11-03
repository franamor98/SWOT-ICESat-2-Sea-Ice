import os
import pandas as pd
import argparse
import xarray as xr
import numpy as np
from tqdm import tqdm

try:
    from shapely.geometry import Polygon, MultiPolygon
    from shapely.ops import unary_union
    _HAS_SHAPELY = True
except Exception:
    _HAS_SHAPELY = False

def _sample_points(lat: np.ndarray, lon: np.ndarray, step_al: int = 10, step_ct: int = 10):
    """
    Return a list of (lon, lat) points from a 2D grid, subsampled both
    along-track and cross-track, ignoring NaNs.
    """
    if lat.ndim != 2 or lon.ndim != 2:
        raise ValueError("Expected 2D lat/lon arrays")
    if lat.shape != lon.shape:
        raise ValueError("lat and lon must have the same shape")

    lat_s = lat[::step_al, ::step_ct]
    lon_s = lon[::step_al, ::step_ct]
    mask = np.isfinite(lat_s) & np.isfinite(lon_s)
    if not np.any(mask):
        return []
    return list(zip(lon_s[mask].ravel().tolist(), lat_s[mask].ravel().tolist()))


def _polygon_wkt_from_points(points: list[tuple[float, float]]):
    """
    Build a single encompassing polygon from points.
    - If shapely is present: convex hull (covers gaps).
    - Else: fallback to bounding box polygon.
    """
    if not points:
        return None

    if _HAS_SHAPELY:
        from shapely.geometry import MultiPoint, Polygon
        hull = MultiPoint(points).convex_hull
        # If hull degenerates to a line/point, use envelope (bbox) to ensure a polygon
        if hull.geom_type == "Polygon":
            return hull.wkt
        else:
            return hull.envelope.wkt
    else:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)
        ring = [(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy), (minx, miny)]
        return f"POLYGON (({', '.join(f'{x} {y}' for x, y in ring)}))"


def get_wkt_swot(file_path: str, step: int = 10):
    """
    Create one polygon that encompasses all valid lat/lon points in the file,
    including gaps between left/right halves or patches.
    Uses convex hull when shapely is available; otherwise falls back to bbox.
    """
    try:
        all_pts = []
        if "Unsmoothed" in file_path and "cropped" not in file_path: #note the raw unsmoothed data its a h5, with left and right groups
            with xr.open_dataset(file_path, group="right") as ds_right:
                all_pts += _sample_points(ds_right["latitude"].values,
                                          ds_right["longitude"].values,
                                          step_al=step, step_ct=step)
            with xr.open_dataset(file_path, group="left") as ds_left:
                all_pts += _sample_points(ds_left["latitude"].values,
                                          ds_left["longitude"].values,
                                          step_al=step, step_ct=step)
        else:
            with xr.open_dataset(file_path) as ds:
                all_pts += _sample_points(ds["latitude"].values,
                                          ds["longitude"].values,
                                          step_al=step, step_ct=step)

        # If subsampling yielded nothing, try without subsampling as a fallback...rarely
        if not all_pts:
            if "Unsmoothed" in file_path and "cropped" not in file_path:
                with xr.open_dataset(file_path, group="right") as ds_right, \
                     xr.open_dataset(file_path, group="left") as ds_left:
                    all_pts += _sample_points(ds_right["latitude"].values,
                                              ds_right["longitude"].values,
                                              step_al=1, step_ct=1)
                    all_pts += _sample_points(ds_left["latitude"].values,
                                              ds_left["longitude"].values,
                                              step_al=1, step_ct=1)
            else:
                with xr.open_dataset(file_path) as ds:
                    all_pts += _sample_points(ds["latitude"].values,
                                              ds["longitude"].values,
                                              step_al=1, step_ct=1)

        return _polygon_wkt_from_points(all_pts)
    except Exception as e:
        print(f"Error building WKT for {file_path}: {e}")
        return None



def parse_swot_filename(swot_file: str, include_wkt: bool, step: int):
    """
    Parse SWOT filename and return its components as a dictionary.
    Expected format:
    SWOT_L2_LR_SSH_Expert_030_412_20250331T234547_20250401T003716_PIC2_01.nc
    """
    try:
        filename = os.path.basename(swot_file)
        parts = filename.split(".")[0].split('_')

        poly_wkt = get_wkt_swot(swot_file, step) if include_wkt else None

        return {
            "product_level": parts[1],        # L2
            "resolution": parts[2],           # LR
            "parameter": parts[3],            # SSH
            "data_type": parts[4],           # Expert
            "cycle": int(parts[5]),          # 030
            "pass": int(parts[6]),           # 412 (ATR)
            "start_time": parts[7],          # 20250331T234547
            "end_time": parts[8],            # 20250401T003716
            "direction": "ascending" if int(parts[6]) % 2 == 1 else "descending",
            "orbit_segment": parts[9],       # PIC2
            "version": parts[10],            # 01
            "polygon_wkt": poly_wkt,         # WKT polygon of the footprint (subsampled)
        }
    except (IndexError, ValueError):
        print(f"Invalid filename format: {filename}")
        return None


def main():
    parser = argparse.ArgumentParser(description="Index SWOT files in a folder with footprint WKT.")
    parser.add_argument("input_folder", type=str, help="Path to the folder containing SWOT .nc files")
    parser.add_argument("--wkt", action="store_true", default=True, help="Include polygon WKT footprint in the output")
    parser.add_argument("--step", type=int, default=10, help="Subsampling step along-track for the footprint (default: 10)")

    args = parser.parse_args()
    input_folder = args.input_folder
    include_wkt = args.wkt
    step = max(1, int(args.step))

    if not os.path.isdir(input_folder):
        print(f"Error: '{input_folder}' is not a valid directory.")
        return

    swot_files = [
        os.path.join(input_folder, f)
        for f in os.listdir(input_folder)
        if f.endswith('.nc')
    ]

    if not swot_files:
        print("No .nc files found in the folder.")
        return

    swot_table = []
    for swot_file in tqdm(swot_files, desc="Processing SWOT files"):
        parsed = parse_swot_filename(swot_file, include_wkt, step)
        if parsed:
            parsed["file_path"] = swot_file
            swot_table.append(parsed)

    if not swot_table:
        print("No valid SWOT files found.")
        return

    df = pd.DataFrame(swot_table)

    # Convert times to datetime
    df["start_time"] = pd.to_datetime(df["start_time"], format="%Y%m%dT%H%M%S")
    df["end_time"] = pd.to_datetime(df["end_time"], format="%Y%m%dT%H%M%S")

    # Sort and reset index
    df.sort_values(by=["cycle", "pass", "start_time"], inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Save to CSV
    output_path = os.path.join(input_folder, "swot_table.csv")
    df.to_csv(output_path, index=False)

    print(f"SWOT table saved to: {output_path}")
    if include_wkt and not _HAS_SHAPELY:
        print("Note: Shapely not found. For 'Unsmoothed' files with left/right groups, the footprint will be a MULTIPOLYGON union (no geometric dissolve). Install shapely for dissolved union.")


if __name__ == "__main__":
    main()
