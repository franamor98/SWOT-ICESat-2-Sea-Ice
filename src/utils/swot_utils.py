import numpy as np
from scipy.interpolate import griddata
import xarray as xr
import pandas as pd
from typing import Literal
from shapely.geometry import LineString, box


def save_field_as_xyz(dataset, field_name, output_base, apply_log=False, clip_range=None):
    """
    Extracts a field from an xarray dataset and saves it as an .xyz file.

    Parameters:
    - dataset (xarray.Dataset): Dataset with 'longitude', 'latitude', and the target field.
    - field_name (str): Name of the field to extract and save.
    - output_base (str): Base name for the output file.
    - apply_log (bool): If True, apply np.log() to the field values.
    - clip_range (tuple or None): If not None, apply np.clip(field, min, max).

    Output:
    - Saves a space-delimited XYZ file with columns: Longitude Latitude Value
    """
    ds = dataset.copy()

    if field_name not in ds.variables:
        raise ValueError(f"Field '{field_name}' not found in dataset.")

    lon = ds['longitude'].values
    lat = ds['latitude'].values
    field = ds[field_name].values

    if apply_log:
        field = np.log(field)

    if clip_range is not None:
        field = np.clip(field, clip_range[0], clip_range[1])

    # Flatten and filter
    lon_flat = lon.flatten()
    lat_flat = lat.flatten()
    field_flat = field.flatten()

    mask = np.isfinite(field_flat)
    df = pd.DataFrame({
        'Longitude': lon_flat[mask],
        'Latitude': lat_flat[mask],
        'Value': field_flat[mask]
    })

    output_file = f"{output_base}"
    df.to_csv(output_file, sep=" ", header=False, index=False)
    print(f"Saved to {output_file}")
    
    return df, output_file


def colocate_swot_is2_tracks(
    swot_table: pd.DataFrame,
    is2_table: pd.DataFrame,
    is2_product = "ATL07",
    swot_time_field: str = "start_time_unsmoothed",
    is2_time_field: str = "date",
) -> pd.DataFrame:
    """
    Time-colocate SWOT and ICESat-2 tracks (ATL03 or ATL07).

    Parameters:
        swot_table (pd.DataFrame): SWOT metadata with at least a time field.
        is2_table (pd.DataFrame): ICESat-2 metadata with a datetime field.
        is2_product (str): Either 'ATL03' or 'ATL07'. Used for suffixing.
        swot_time_field (str): SWOT column to use for time comparison.
        is2_time_field (str): ICESat-2 column to use for time comparison.
        top_n (int): Number of closest matches to return.

    Returns:
        pd.DataFrame: Sorted colocated matches with time difference.
    """

    if is2_product not in {"ATL03", "ATL07"}:
        raise ValueError("is2_product must be either 'ATL03' or 'ATL07'")

    # Ensure datetimes
    swot_table = swot_table.copy()
    is2_table = is2_table.copy()
    swot_table[swot_time_field] = pd.to_datetime(swot_table[swot_time_field])
    is2_table[is2_time_field] = pd.to_datetime(is2_table[is2_time_field])

    swot_table.sort_values(by=swot_time_field, inplace=True)
    is2_table.sort_values(by=is2_time_field, inplace=True)
    is2_table = is2_table[is2_table[is2_time_field].notna()]
    is2_table = is2_table[is2_table["init_lat"].notna()]

    swot_table.reset_index(drop=True, inplace=True)
    is2_table.reset_index(drop=True, inplace=True)

    # Merge on time using nearest strategy
    merged = pd.merge_asof(
        swot_table,
        is2_table,
        left_on=swot_time_field,
        right_on=is2_time_field,
        direction='nearest',
        suffixes=('_swot', f'_{is2_product.lower()}')
    )

    # Compute time difference
    merged['time_diff'] = (merged[is2_time_field] - merged[swot_time_field]).abs()
    
    #rename filepath to filepath_is2
    merged.rename(columns={"filepath": f"file_path_{is2_product.lower()}"}, inplace=True)

    
    #determine which pairs are space located 
    merged["is2_segment"] = merged.apply(
        lambda row: LineString([
            (row['init_lon'], row['init_lat']),
            (row['end_lon'], row['end_lat'])
        ]) 
    , axis=1)

    merged["swot_bbox"] = merged.apply(
        lambda row: box(
            row['min_lon'], row['min_lat'],
            row['max_lon'], row['max_lat']
        ), axis=1
    )

    # Check if the ICESat-2 segment intersects with the SWOT bounding box
    merged["intersects"] = merged.apply(
        lambda row: row["is2_segment"].intersects(row["swot_bbox"]), axis=1
    )

    # drop segment and bbox columns
    merged.drop(columns=["is2_segment", "swot_bbox"], inplace=True)
    # get only columns of intersecting pairs
    merged = merged[merged["intersects"]]
    
  
    # Return top N matches sorted by time difference
    return merged.sort_values('time_diff').reset_index(drop=True)





def crop_swot_dataset(ds, lat_bounds=(70, 79), lon_bounds=(-150, -130)):
    """
    Convert longitude from 0–360 to -180–180, then crop to lat/lon bounds.
    If no data remains after cropping, return None.
    """
    # Convert longitude
    ds['longitude'] = ds['longitude'].where(ds['longitude'] <= 180, ds['longitude'] - 360)

    # Unpack bounds
    lat_min, lat_max = lat_bounds
    lon_min, lon_max = lon_bounds

    # Mask
    mask = (
        (ds.latitude >= lat_min) & (ds.latitude <= lat_max) &
        (ds.longitude >= lon_min) & (ds.longitude <= lon_max)
    )
    ds_cropped = ds.where(mask, drop=True)

    # Check if anything is left (e.g. num_lines > 0)
    if ds_cropped.sizes.get("num_lines", 0) == 0:
        return None

    return ds_cropped



def load_unsmoothed_swot(file_path, lat_bounds=(70, 79), lon_bounds=(-150, -130)):
    """
    Load left and right SWOT data, merge them, and crop to the specified bounds.
    """
    swot_left = xr.open_dataset(file_path, group="left")
    swot_right = xr.open_dataset(file_path, group="right")
    
    # Sanity check
    assert set(swot_left.data_vars) == set(swot_right.data_vars), "Variables differ between left and right"
    
    # Concatenate
    swot_combined = xr.concat([swot_left, swot_right], dim="num_lines")

    
    # Crop
    swot_cropped = crop_swot_dataset(swot_combined, lat_bounds, lon_bounds)
    
    swot_left.close()
    swot_right.close()
    
    return swot_cropped



def interpolate_unsmoothed(
        expert_ds: xr.Dataset, unsmoothed_ds: xr.Dataset, 
        vars_to_interpolate: list
        ) -> xr.Dataset:
    """
    Interpolate selected variables from expert_ds onto the unsmoothed_ds grid.

    Parameters:
    - expert_ds (xr.Dataset): Source dataset with expert corrections.
    - unsmoothed_ds (xr.Dataset): Target dataset where variables will be interpolated.
    - vars_to_interpolate (list): List of variable names to interpolate.

    Returns:
    - xr.Dataset: Updated unsmoothed_ds with interpolated variables.
    """
    lon_src = expert_ds["longitude"].values.flatten()
    lat_src = expert_ds["latitude"].values.flatten()
    points_src = np.stack([lon_src, lat_src], axis=-1)

    lon_tgt = unsmoothed_ds["longitude"].values
    lat_tgt = unsmoothed_ds["latitude"].values
    
    expert_ds["rad_wet_tropo_cor"] = expert_ds["rad_wet_tropo_cor"].fillna(0)
    expert_ds["wet_tropo_cor"] = expert_ds.rad_wet_tropo_cor-expert_ds.model_wet_tropo_cor
    vars_to_interpolate.append("wet_tropo_cor")

    for var in vars_to_interpolate:
        print(f"Interpolating '{var}' → unsmoothed grid")
        values_src = expert_ds[var].values.flatten()

        values_interp = griddata(
            points_src, values_src,
            (lon_tgt, lat_tgt),
            method="linear",
            fill_value=np.nan
        )

        unsmoothed_ds[var + "_expert"] = (("num_lines", "num_pixels"), values_interp)

    return unsmoothed_ds







def apply_ssha_corrections(
    ds: xr.Dataset,
    mss_field: str = "mean_sea_surface_cnescls",
    return_expert: bool = False
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Compute SSHA high-resolution (HR) and optionally low-resolution (LR) corrections.
    Parameters:
    - ds (xr.Dataset): Dataset with required fields.
    - mss_field (str): Base name of the mean sea surface field (default: 'mean_sea_surface_cnescls').
    - return_expert (bool): Whether to compute and return the LR (expert-interpolated) SSHA.

    Returns:
    - tuple:
        - ssha_HR (np.ndarray): High-resolution SSHA array.
        - ssha_LR (np.ndarray or None): Low-resolution SSHA array (if return_expert=True), else None.
    """

    
    wet_tropo_cor = ds["wet_tropo_cor"] if "wet_tropo_cor" in ds else 0
    
    ssha_HR = (
        ds["ssh_karin_2"]
        + wet_tropo_cor
        + ds["height_cor_xover_expert"]
        - ds[mss_field]
        - ds["solid_earth_tide_expert"]
        - ds["ocean_tide_fes_expert"]
        - ds["ocean_tide_non_eq_expert"]
        - ds["internal_tide_hret_expert"]
        - ds["pole_tide_expert"]
        - ds["dac_expert"]
    )

    ssha_LR = None
    if return_expert:
        ssha_LR = (
            ds["ssh_karin_2_expert"]
            + wet_tropo_cor
            + ds["height_cor_xover_expert"]
            - ds[mss_field]
            - ds["solid_earth_tide_expert"]
            - ds["ocean_tide_fes_expert"]
            - ds["ocean_tide_non_eq_expert"]
            - ds["internal_tide_hret_expert"]
            - ds["pole_tide_expert"]
            - ds["dac_expert"]
        )

    return ssha_HR, ssha_LR