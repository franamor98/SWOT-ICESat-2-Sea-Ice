
import argparse
import os
from pathlib import Path
import xarray as xr
from tqdm import tqdm

import sys
sys.path.append(".")
from src.utils.swot_utils import crop_swot_dataset, load_unsmoothed_swot

def main():
    parser = argparse.ArgumentParser(description="Process SWOT .nc files with optional cropping.")
    parser.add_argument("input_dir", type=str, help="Path to the folder containing .nc files")
    parser.add_argument("--lat_bounds", nargs=2, type=float, default=[70, 79],
                        help="Latitude bounds (default: 70 79)")
    parser.add_argument("--lon_bounds", nargs=2, type=float, default=[-150, -130],
                        help="Longitude bounds (default: -150 -130)")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = input_dir.parent / f"{input_dir.name}_cropped"
    output_dir.mkdir(exist_ok=True)
    print(f"Output directory: {output_dir}")

    lat_bounds = tuple(args.lat_bounds)
    lon_bounds = tuple(args.lon_bounds)

    is_unsmoothed = "unsmoothed" in input_dir.name.lower()

    for file_path in tqdm(input_dir.glob("*.nc"), desc="Processing files"):

        try:

            bbox_str = f"{lat_bounds[0]}_{lat_bounds[1]}_{lon_bounds[0]}_{lon_bounds[1]}"
            new_filename = file_path.stem + bbox_str + file_path.suffix
            output_path = output_dir / new_filename

            if output_path.exists():
                tqdm.write(f"\tSkipping {file_path.name} because {output_path.name} already exists.")
                continue

            if is_unsmoothed:
                ds = load_unsmoothed_swot(file_path, lat_bounds, lon_bounds)
            else:
                with xr.open_dataset(file_path) as dataset:
                    ds = crop_swot_dataset(dataset, lat_bounds, lon_bounds)


            ds.to_netcdf(output_path)
            ds.close()
        except Exception as e:
            tqdm.write(f"\tError processing {file_path}: {e}")
            continue

if __name__ == "__main__":
    main()