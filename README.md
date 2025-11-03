# README — SWOT–ICESat-2 Matches
Author: Fran Amor  
Affiliation: UiT The Arctic University of Norway, Tromsø  

## Overview
This repository provides tools and workflows to analyze matching passes of ICESat-2 (ATL03) and SWOT (Surface Water and Ocean Topography) over sea ice.

It contains data download, subsetting, and co-location of SWOT and ICESat-2 measurements to support sea-ice elevation and SSH anomaly analyses.

An example web visualizatin of matching measurements can be found here:
https://franamor98.github.io/Sea_ice_Viz/#8/73.053/-139.614

## Repository Structure
```bash

.
├── config
│   ├── Is2_conf.json
│   └── swot_conf.json
├── data
├── environment.yml
├── notebooks
│   └── swot_IS2.ipynb
└── src
    ├── data_services
    │   ├── Icesat2_data_service_ea.py
    │   └── swot_data_service.py
    ├── is2_table.py
    ├── swot_table.py
    └── utils
        ├── Is2_utils.py
        ├── process_swot_files.py
        └── swot_utils.py
```

## Installation
It is recommended to use conda:

```bash
# Clone repository
git clone https://github.com/youruser/swot_is2_matches.git
cd swot_is2_matches

# Create environment
conda env create -f environment.yml
conda activate swot-is2
```

## Configuration
You define your download parameters in the JSON configuration files inside the config folder.

### Example: ICESat-2 configuration (Is2_conf.json)
```json
{
  "product": "ATL03",
  "spatial_extent": [-150, 70, -130, 76],
  "date_range": ["2024-02-01", "2024-04-16"],
  "download_path": "./data/IS2/ATL03_2"
}
```

### Example: SWOT configuration (swot_conf.json)
```json
{
  "product": "SWOT_L2_LR_SSH_2.0",
  "spatial_extent": [-150, 70, -130, 76],
  "date_range": ["2024-02-01", "2024-02-15"],
  "filter": "Unsmoothed",
  "download_path": "./data/SWOT/L2_LR_SSH_2.0/"
}
```

Note:
- `filter` can be "Unsmoothed" or "Expert", depending on which SWOT L2 LR sub-product you are downloading.
- It is recommended to download both types in separate folders:
  - Unsmoothed: main geophysical measurements.
  - Expert: auxiliary corrections (tidal, atmospheric, mean sea surface).
- The Expert product is lighter and can be downloaded in batch alongside the Unsmoothed data.

## Processing Workflow

### 1. Crop SWOT passes to your AOI
The SWOT data service stores full-orbit passes. You can subset them spatially using the process_swot_files.py utility:

```bash
python src/utils/process_swot_files.py ./data/SWOT/L2_LR_SSH_2.0/Unsmoothed     --lat_bounds 70 79     --lon_bounds -150 -130
```

This creates cropped SWOT files for your Area of Interest (AOI). Repeat this step for both Unsmoothed and Expert folders.

### 2. Build index tables
Generate index tables (metadata + geometry polygons) for each dataset:

```bash
python src/swot_table.py ./data/SWOT/L2_LR_SSH_2.0/Unsmoothed_cropped
python src/swot_table.py ./data/SWOT/L2_LR_SSH_2.0/Expert_cropped
python src/is2_table.py ./data/IS2/ATL03_2
```

Each script will output a *_table.csv file containing temporal and spatial metadata and a simplified WKT polygon representing the pass footprint.

### 3. Match SWOT and ICESat-2 passes
Use the notebook notebooks/swot_IS2.ipynb to:
- Load the index CSVs.
- Match Unsmoothed–Expert SWOT pairs.
- Match SWOT passes with ICESat-2 granules using time and geometry.

Temporal matching logic:
- SWOT orbit segments take roughly ~50 minutes per pass. File provides end and start time of each half orbit (pass)
- For high-latitude regions:
  - Ascending passes: Arctic data correspond to the end time of the SWOT file.
  - Descending passes: Arctic data correspond to the start time.
- ICESat-2 ATL03 granules last about 7 minutes (1/14 of an orbit)

Overall... lets say time uncertainty is approximate within ±10 minutes. For the sake of simplicity this non robust method is kept.

### 4. Colocation and analysis
The notebook then:
1. Loads both Unsmoothed and Expert SWOT data.
2. Interpolates geophysical corrections from Expert to Unsmoothed.
3. Applies corrections to obtain SSHA (Sea Surface Height Anomaly).
4. Loads ATL03 photons, applies geophysical corrections, and smooths them (500-sample moving window, using photons with sea_ice_conf > 1).
5. Finds ICESat-2 photons within the SWOT swath and computes nearest-neighbor SWOT–IS2 matches.

Reference for processing approach:
Reint Fischer, Sinead L. Farrell, Kyle Duncan, et al. (2025).
Swath Mapping Altimetry over Sea Ice: First Results from the Surface Water and Ocean Topography (SWOT) Mission.
ESS Open Archive, May 2025.



Author: Fran Amor  
UiT – The Arctic University of Norway, Tromsø