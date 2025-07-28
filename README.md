# teresa
Code supporting MicroRider microstructure data processing
# Glider Microstructure Turbulence Processing

This repository contains MATLAB and Python scripts used for processing microstructure data collected by the RSI MicroRider mounted on the Slocum glider "Teresa" (SMART missions, 2015–2024, Western Mediterranean). The data include turbulence dissipation rate estimates (ε, χ) from shear and thermistor sensors.

## Description

The code implements best practices for ocean turbulence processing as described by:
- Lueck et al. (2024) [https://doi.org/10.3389/fmars.2024.1334327](https://doi.org/10.3389/fmars.2024.1334327)
- Piccolroaz et al. (2021)
- ODAS Matlab Toolbox v4.5.1 (Rockland Scientific)

It includes:
- Preprocessing of L0-L1 glider microstructure files
- Spectral integration to estimate ε and χ
- Quality Control (QC) flags
- Section and profile-based netCDF export
- Optional integration with flight model by Merckelbach et al. (2019)

## License

This code is released under the [MIT License](LICENSE), allowing unrestricted use, modification, and redistribution.

## How to cite

If you use this code, please cite:

Kokoszka, F. V. M., et al. (2025). Glider Microstructure Turbulence Processing Code. Zenodo. 10.5281/zenodo.16541936
https://zenodo.org/records/16541936
