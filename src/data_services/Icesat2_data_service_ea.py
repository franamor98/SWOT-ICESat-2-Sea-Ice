import earthaccess as ea
import json
import argparse
import os
from datetime import datetime


def load_config(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config


def query_and_download(config):
    product = config.get("product")
    if not product:
        raise ValueError("Product must be specified in the config file.")

    spatial_extent = config.get("spatial_extent", [-180, -90, 180, 90])
    date_range = config.get("date_range", ["2023-01-01", "2023-01-31"])
    download_path = config.get("download_path", f"./IS2/{product}")
    os.makedirs(download_path, exist_ok=True)

    # Unpack bounding box
    if len(spatial_extent) != 4:
        raise ValueError("Spatial extent must be a list of four values: [min_lon, min_lat, max_lon, max_lat]")

    min_lon, min_lat, max_lon, max_lat = spatial_extent
    start_date, end_date = date_range

    print(f"Querying {product} from {start_date} to {end_date} in region {spatial_extent}...")

    results = ea.search_data(
        short_name=product,
        temporal=(start_date, end_date),
        bounding_box=(min_lon, min_lat, max_lon, max_lat),
        cloud_hosted=False
    )

    if not results:
        print("No granules found.")
        return

    print(f"Found {len(results)} granules. Downloading...")

    downloaded = ea.download(results, download_path)

    print(f"Downloaded {len(downloaded)} files to {download_path}.")


def main():
    parser = argparse.ArgumentParser(description="Download ICESat-2 data using earthaccess")
    parser.add_argument("--config", default="./config/Is2_conf.json", help="Path to JSON config file")

    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Config file {args.config} not found.")
        return

    print("Authenticating with Earthdata...")
    ea.login()

    print(f"Loading configuration from {args.config}...")
    config = load_config(args.config)
    query_and_download(config)


if __name__ == "__main__":
    main()
