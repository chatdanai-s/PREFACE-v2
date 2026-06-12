# PREFACE
An Automated Pipeline for the Selection of Transmission Spectroscopy Candidates

This repository is work in process, here is the to-do list:
* [ ] Refactor PREFACE in pip-installable form (Phase 1 complete)
* [ ] Increase code maintainability by giving human-readable flag names
* [x] Include Kempton et al. (2018)'s TSM score in the output
* [x] Update Mass-Radius(-Temperature) relation according to Edmondson et al. (2023)
* Remove LIME_tbx and AERONET dependence for moonlight scattering modeling
  - Completely use lookup tables for lunar altaz and magnitude time series. The lookup tables will be in 2 year intervals for 2026-2050, and they must be stored in float32 and compressed for space conservation.
  - Leave LUT and aggregated AERONET generation code in case of future LUT updates
  - LIME_tbx still referenced for empirical lunar magnintude conversion between filters
