import earthaccess as ea
import json
import argparse
import os


def load_config(config_file):
    with open(config_file, 'r') as f:
        config = json.load(f)
    return config

def filter_granules(granules, filter_type):

    return [g for g in granules if filter_type in g['meta']['native-id']]
    


def query_and_download(config):
    product = config.get("product", "SWOT_L2_LR_SSH_2.0")
    spatial_extent = config.get("spatial_extent", [-150, 70, -130, 79])
    date_range = config.get("date_range", ["2025-04-01", "2025-04-30"])
    filter_type = config.get("filter", "")
    download_path = config.get("download_path", f"./SWOT/{product}_{filter_type.lower()}")

    os.makedirs(download_path, exist_ok=True)

    if len(spatial_extent) != 4:
        raise ValueError("Spatial extent must be [min_lon, min_lat, max_lon, max_lat]")

    min_lon, min_lat, max_lon, max_lat = spatial_extent
    start_date, end_date = date_range

    print(f"Querying {product} ({filter_type}) from {start_date} to {end_date} in {spatial_extent}...")

    results = ea.search_data(
        short_name=product,
        temporal=(start_date, end_date),
        #bounding_box=(min_lon, min_lat, max_lon, max_lat),
    )
    print(f"Found {len(results)} granules before filtering.")    
    filtered_results = filter_granules(results, filter_type)
    print(f"Found {len(filtered_results)} granules after filtering with '{filter_type}'.")

    if not filtered_results:
        print("No granules found with the specified filter.")
        return

    print(f"Found {len(filtered_results)} granules. Downloading to {download_path}...")

    downloaded = ea.download(filtered_results, download_path)
    print(f"Downloaded {len(downloaded)} files.")


def main():
    parser = argparse.ArgumentParser(description="Download SWOT data using earthaccess")
    parser.add_argument("--config", default="./config/swot_conf.json", help="Path to JSON config file")
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
