import h5py
import pandas as pd
import numpy as np
from typing import Tuple


def apply_moving_average(df, column='height', window=5, group_by='beam'):
    """
    Returns a copy of the DataFrame with a moving average applied to the specified column,
    grouped by the specified field (default 'beam').

    Parameters:
        df (pd.DataFrame): Input DataFrame with at least 'beam' and the target column.
        column (str): The column to apply the moving average to.
        window (int): Window size for the moving average.
        group_by (str): The column to group by (default 'beam').

    Returns:
        pd.DataFrame: Copy of the input DataFrame with a new column `{column}_smoothed`.
    """
    df_copy = df.copy()
    df_copy[f"{column}_smoothed"] = (
        df_copy
        .groupby(group_by)[column]
        .transform(lambda x: x.rolling(window, center=True, min_periods=1).mean())
    )
    return df_copy



def find_colocated_ATL_files(
    file: str,
    table_07: pd.DataFrame,
    table_03: pd.DataFrame
) -> Tuple[str, str, pd.Timedelta]:
    """
    Given a filename from ATL07 or ATL03, return the closest matching file
    from the other product based on date proximity.

    Parameters:
        file (str): Filename of ATL07 or ATL03 product.
        table_07 (pd.DataFrame): Table with ATL07 metadata (must include 'filename' and 'date').
        table_03 (pd.DataFrame): Table with ATL03 metadata (must include 'filename' and 'date').

    Returns:
        Tuple[str, str, pd.Timedelta]: (ATL07 filename, ATL03 filename, time difference)
    """
    # Ensure datetime
    table_07 = table_07.copy()
    table_03 = table_03.copy()
    table_07["date"] = pd.to_datetime(table_07["date"])
    table_03["date"] = pd.to_datetime(table_03["date"])
    table_07["day"] = table_07["date"].dt.date
    table_03["day"] = table_03["date"].dt.date

    # Identify product type
    if "ATL07" in file:
        reference_table = table_07
        other_table = table_03
    elif "ATL03" in file:
        reference_table = table_03
        other_table = table_07
        
    else:
        raise ValueError("Filename must contain 'ATL07' or 'ATL03'")

    # Extract date for the input file
    ref_row = reference_table[reference_table["filepath"] == file]
    if ref_row.empty:
        raise ValueError(f"File {file} not found in its respective table.")

    ref_date = ref_row.iloc[0]["date"]
    ref_day = ref_date.date()

    # Filter other table by day
    same_day_matches = other_table[other_table["day"] == ref_day].copy()
    if same_day_matches.empty:
        raise ValueError(f"No same-day matches found for {file}")

    same_day_matches["time_diff"] = (same_day_matches["date"] - ref_date).abs()
    best_match = same_day_matches.sort_values("time_diff").iloc[0]

    if "ATL07" in file:
        return file, best_match["filepath"], best_match["time_diff"]
    else:
        return best_match["filepath"], file, best_match["time_diff"]




def get_beam_strengths(file_path):
    """
    Determine the strong and weak beams in an ICESat-2 HDF5 file based on spacecraft orientation.

    Parameters:
    file_path (str): Path to the ICESat-2 HDF5 file.

    Returns:
    tuple: Lists of strong and weak beams.
    """

    with h5py.File(file_path, 'r') as f:
        # Get orientation: 0 = backward (left beams strong), 1 = forward (right beams strong)
        sc_orient = f["/orbit_info/sc_orient"][0]

        beam_pairs = {
            "gt1l": "Left", "gt1r": "Right",
            "gt2l": "Left", "gt2r": "Right",
            "gt3l": "Left", "gt3r": "Right"
        }

        strong_beams = []
        weak_beams = []

        for beam, side in beam_pairs.items():
            if beam in f:
                is_strong = (sc_orient == 0 and side == "Left") or (sc_orient == 1 and side == "Right")
                if is_strong:
                    strong_beams.append(beam)
                else:
                    weak_beams.append(beam)

      
        f.close()

        return strong_beams, weak_beams
    

import numpy as np
import pandas as pd
import h5py

def load_ATL03_file_with_interp(
    filepath,
    date,
    track_id=None,
    min_conf_ocean=4,
    min_conf_seaice=4,
    max_conf_land=0,
):
    """
    Loads ATL03 photon height data from a single file, applies geophysical interpolation,
    and returns a memory-efficient DataFrame.

    Parameters:
        filepath (str): Path to ATL03 file.
        date (str): Date string to tag data.
        track_id (str or int, optional): Track ID label.
        min_conf_ocean (int): Minimum signal confidence for ocean photons (0–4).
        min_conf_seaice (int): Minimum signal confidence for sea ice photons (0–4).
        max_conf_land (int): Maximum allowable land signal confidence (0–4).
    """
    try:
        strong_beams, weak_beams = get_beam_strengths(filepath)
        beams = strong_beams + weak_beams
    except Exception as e:
        print(f"Skipping {filepath} due to beam detection error: {e}")
        return pd.DataFrame()

    beam_dfs = []

    try:
        with h5py.File(filepath, "r") as f:
            for beam in beams:
                try:
                    beam_type = "strong" if beam in strong_beams else "weak"
                    group = f[beam]
                    heights = group["heights"]
                    geophys = group["geophys_corr"]
                    geoloc = group["geolocation"]

                    # Load variables
                    lat = heights["lat_ph"][:]
                    lon = heights["lon_ph"][:]
                    h_ph = heights["h_ph"][:]
                    delta_time = heights["delta_time"][:]
                    quality = heights["quality_ph"][:]
                    signal_conf = heights["signal_conf_ph"][:]  # shape (n, 4)

                    # Extract relevant columns
                    signal_conf_land = signal_conf[:, 0]
                    signal_conf_ocean = signal_conf[:, 1]
                    signal_conf_seaice = signal_conf[:, 2]

                    # Apply confidence mask if needed
                    keep_mask = (
                        ((signal_conf_seaice >= min_conf_seaice) |
                         (signal_conf_ocean >= min_conf_ocean)) &
                        (signal_conf_land <= max_conf_land)
                    )

                    if not np.any(keep_mask):
                        continue

                    # Apply mask and cast types
                    lat = lat[keep_mask].astype(np.float32)
                    lon = lon[keep_mask].astype(np.float32)
                    h_ph = h_ph[keep_mask].astype(np.float32)
                    delta_time = delta_time[keep_mask]
                    quality = quality[keep_mask].astype(np.int8)
                    signal_conf_land = signal_conf_land[keep_mask].astype(np.int8)
                    signal_conf_ocean = signal_conf_ocean[keep_mask].astype(np.int8)
                    signal_conf_seaice = signal_conf_seaice[keep_mask].astype(np.int8)

                    # Interpolate corrections (on filtered delta_time)
                    ref_dt = geoloc["delta_time"][:]
                    distance_x = geoloc["segment_dist_x"][:]
                    geoid = np.interp(delta_time, ref_dt, geophys["geoid"][:]).astype(np.float32)
                    dac = np.interp(delta_time, ref_dt, geophys["dac"][:]).astype(np.float32)
                    tide_ocean = np.interp(delta_time, ref_dt, geophys["tide_ocean"][:]).astype(np.float32)
                    tide_load = np.interp(delta_time, ref_dt, geophys["tide_load"][:]).astype(np.float32)
                    distance_x = np.interp(delta_time, ref_dt, distance_x).astype(np.float32)

                    n = len(lat)

                    df = pd.DataFrame({
                        "latitude": lat,
                        "longitude": lon,
                        "height": h_ph,
                        "delta_time": delta_time,
                        "signal_conf_land": signal_conf_land,
                        "signal_conf_ocean": signal_conf_ocean,
                        "signal_conf_seaice": signal_conf_seaice,
                        "filepath": np.full(n, filepath),
                        "quality": quality,
                        "beam": np.full(n, beam),
                        "beam_type": np.full(n, beam_type),
                        "geoid": geoid,
                        "dac": dac,
                        "tide_ocean": tide_ocean,
                        "load_tide": tide_load,
                        "date": np.full(n, date),
                        "distance_x": distance_x,
                    })

                    if track_id is not None:
                        df["track_id"] = np.full(n, track_id)

                    beam_dfs.append(df)

                except KeyError as e:
                    print(f"Missing variable in beam {beam} of {filepath}: {e}")
                    continue

    except OSError as e:
        print(f"Could not open file {filepath}: {e}")

    return pd.concat(beam_dfs, axis=0, ignore_index=True, copy=False) if beam_dfs else pd.DataFrame()



def load_ATL07_file(filepath, date, track_id=None):
    """
    Loads sea ice segment data from a single ATL07 file.

    Parameters:
        filepath (str): Path to the ATL07 HDF5 file.
        date (datetime): Timestamp associated with the file.
        track_id (str or None): Optional track ID.

    Returns:
        pd.DataFrame: DataFrame containing sea ice segment data for all beams.
    """
    try:
        strong_beams, weak_beams = get_beam_strengths(filepath)
        beams = strong_beams + weak_beams
    except Exception as e:
        print(f"Skipping {filepath} due to beam detection error: {e}")
        return pd.DataFrame()

    beam_dfs = []

    try:
        with h5py.File(filepath, "r") as f:
            for beam in beams:
                try:
                    beam_type = "strong" if beam in strong_beams else "weak"

                    lat = f[beam]["sea_ice_segments"]["latitude"][:]
                    lon = f[beam]["sea_ice_segments"]["longitude"][:]
                    delta_t = f[beam]["sea_ice_segments"]["delta_time"][:]
                    height = f[beam]["sea_ice_segments"]["heights"]["height_segment_height"][:]
                    ph_type = f[beam]["sea_ice_segments"]["heights"]["height_segment_type"][:]
                    n_pulse_seg = f[beam]["sea_ice_segments"]["heights"]["height_segment_n_pulse_seg"][:]
                    quality = f[beam]["sea_ice_segments"]["heights"]["height_segment_quality"][:]
                    confidence = f[beam]["sea_ice_segments"]["heights"]["height_segment_confidence"][:]
                    distance = f[beam]["sea_ice_segments"]["seg_dist_x"][:]
                    fit_quality_flag = f[beam]["sea_ice_segments"]["heights"]["height_segment_fit_quality_flag"][:]

                    df = pd.DataFrame({
                        "latitude": lat,
                        "longitude": lon,
                        "height": height,
                        "ph_type": ph_type,
                        "n_pulse_seg": n_pulse_seg,
                        "quality": quality,
                        "confidence": confidence,
                        "beam": beam,
                        "beam_type": beam_type,
                        "date": date,
                        "distance": distance,
                        "fit_quality_flag": fit_quality_flag,
                        "delta_time": delta_t
                    })

                    if track_id is not None:
                        df["track_id"] = track_id

                    beam_dfs.append(df)

                except KeyError as e:
                    print(f"Missing variable in beam {beam} of {filepath}: {e}")
                    continue

    except OSError as e:
        print(f"Could not open file {filepath}: {e}")

    return pd.concat(beam_dfs, ignore_index=True) if beam_dfs else pd.DataFrame()




def load_ATL07_tracks(
    subset: pd.DataFrame,
    filepath_field: str = "filepath",
    track_id_field: str = "track_id"
) -> pd.DataFrame:
    """
    Loads ATL07 data for each row in the DataFrame using `load_ATL07_file`.

    Returns:
        pd.DataFrame: Concatenated data from all files.
    """
    all_data = []

    for _, row in subset.iterrows():
        filepath = row[filepath_field]
        date = pd.to_datetime(row["date"])
        track_id = row[track_id_field] if track_id_field and track_id_field in row else None

        df = load_ATL07_file(filepath, date, track_id)
        if not df.empty:
            all_data.append(df)

    return pd.concat(all_data, ignore_index=True) if all_data else pd.DataFrame()

