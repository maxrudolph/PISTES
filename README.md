# PISTES - Planetary Ice Shell Thermal Evolution and Stress
## Thermal and stress evolution of icy moons, with application to Europa, Enceladus, Mimas, and other Icy Ocean Worlds.


### Papers that use PISTES:
1. Rudolph, M.L., Manga, M., Walker, M., and Rhoden, A.R. (2022). Cooling crusts create concomittant cryovolcanic cracks. Geophysical Research Letters 49(5),  e2021GL094421.
2. Rhoden, A. R., Rudolph, M. L., & Manga, M. (2023). The challenges of driving Charon's cryovolcanism from a freezing ocean. Icarus, 392, 115391.
3. Rhoden, A.R., Walker, M.E., Rudolph, M.L., Bland, M.T., and Manga, M. (2024) The evolution of a young ocean within Mimas. Earth and Planetary Science Letters 635, 118698.
4. Rudolph, M.L., Manga, M., Rhoden, A.R., and Walker, M. (2025) Boiling oceans and compressional tectonics on emerging ocean worlds. Nature Astronomy.

### Code citation:
If you use PISTES in your own work, please cite Rudolph et al. (2022), Rhoden et al. (2024), and Rudolph et al. (2025).

### Requirements:
- MATLAB 2020a or later due to the use of tiled layouts.
- Fabio Crameri's scientific colormaps (https://www.fabiocrameri.ch/colourmaps/)
- Parallel computing toolbox for parfor loop parallelization (optional)

### Description of files:
- The equations solved by the code are derived in the Supplemental Material to Rudolph et al. (2022): ```2021GL094421_Supplemental_Material.pdf```
- The main driver program that sweeps parameter space to reproduce the figures in Rudolph et al. (2022) is ```sweep_parameter_space.m```. This will produce about 16GB of output as .mat files.
- The file ```thickening_ice_shell.m``` solves for the stresses in an initially stress-free thickening ice shell with no heat input.
- Helper subroutines are located in the subdirectory ```core/```
- Files in the ```benchmarks/``` subdirectory can be used to reproduce the results in Nimmo (2004).
- Files in the ```underpressure/``` folder can be used to reproduce the results of Rudolph et al. (2025)
    - ```Rudolph_etal_2025_Figure3_underpressure_regime_plot.m``` - reproduce Figure 3 in the main text and Figures S1-S5 in the Supplementary Information
    - ```Rudolph_etal_2025_Figure4_enceladus_mimas_thinning_underpressure.m``` - reproduce Figure 4 in the main text and Figures S6-7 in the Supplementary Information

