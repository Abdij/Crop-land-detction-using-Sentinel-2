# Crop-land-detction-using-Sentinel-2
# Farmland Detection in Somalia using Sentinel-2

A Google Earth Engine (GEE) script for detecting and visualizing agricultural land in Somalia using Sentinel-2 satellite imagery and vegetation indices.

## Overview

This project identifies likely farmland areas across major agricultural regions in Somalia using multi-index analysis from Sentinel-2 imagery. It combines vegetation, moisture, and soil-related indices to classify cultivated land and filter out non-agricultural areas such as bare soil, urban areas, water, and dense forest.

The script is designed for agricultural mapping, land monitoring, and spatial analysis in arid and semi-arid environments.

## Study Regions

The script includes predefined areas for major farming zones in Somalia:

- Lower Shabelle
- Middle Shabelle
- Bay Region
- Juba Valley
- Gabiley

You can switch the analysis area by changing the `studyArea` variable in the script.

## Data Source

- **Satellite Data:** Sentinel-2 Surface Reflectance Harmonized
- **Dataset ID:** `COPERNICUS/S2_SR_HARMONIZED`
- **Temporal Coverage:** 2024, with fallback to 2023 if no suitable scenes are found
- **Spatial Resolution:** 10 meters

## Methodology

The script performs the following steps:

1. Defines the study area
2. Loads Sentinel-2 imagery for the selected region
3. Applies cloud masking using the `QA60` band
4. Builds a median composite image
5. Calculates multiple spectral indices:
   - NDVI (Normalized Difference Vegetation Index)
   - EVI (Enhanced Vegetation Index)
   - SAVI (Soil-Adjusted Vegetation Index)
   - NDMI (Normalized Difference Moisture Index)
   - BSI (Bare Soil Index)
6. Applies farmland classification rules
7. Removes isolated small patches using connected pixel filtering
8. Displays results on the map with a custom legend
9. Calculates farmland statistics
10. Exports farmland and NDVI layers to Google Drive

## Farmland Detection Criteria

A pixel is classified as farmland when it satisfies all of the following conditions:

- `NDVI > 0.2` → vegetated area
- `NDVI < 0.7` → excludes dense forest
- `SAVI > 0.15` → suitable vegetation in arid environments
- `NDMI > -0.1` → some moisture present
- `BSI < 0` → excludes bare soil
- Minimum patch size of **0.5 hectares**

## Visualization Layers

The script adds several map layers for interpretation:

- Study Area boundary
- Natural Color Composite (RGB: B4, B3, B2)
- False Color Composite (B8, B4, B3)
- Agriculture Composite (B11, B8, B4)
- NDVI
- SAVI
- NDMI
- Final detected farmland areas

## Outputs

The script produces:

- A farmland detection layer
- Vegetation and moisture index layers
- Area statistics in hectares
- Percentage of farmland in the selected study area
- Export tasks for:
  - Farmland raster
  - NDVI raster

## Requirements

To use this project, you need:

- A [Google Earth Engine](https://earthengine.google.com/) account
- Access to the GEE Code Editor

## How to Use

1. Open the Google Earth Engine Code Editor
2. Copy and paste the script into a new script file
3. Select your preferred study region by modifying:

```javascript
var studyArea = lowerShabelle;
