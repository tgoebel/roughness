# roughness

Compute roughness based on 2D-surface scan. The contained example file is a white-light interferometry scan
of micron-scale lab-fault topography.

Analysis steps:
1) create regular grid and remove outliers, fill data gaps
2) detrend the surface
3) stack slip parallel or orthogomal profiles
4) fit the PSD of stacked profiles and compute exponent based on lsq of log-transformed data


Input file: fault_roughness.txt
!make sure to adjust dir_in to point to right data directory before running the script
