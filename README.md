# logKcalc

This R package calculates equilibrium constants for thermodynamic data files used in [The Geochemist's Workbench®](https://www.gwb.com).
It reads reactions from an input GWB data file, calculates equilibrium constants (log&#8201;*K*) for the reactions using the [OBIGT database in CHNOSZ](http://chnosz.net/vignettes/obigt.html), and writes them to an output GWB file.

## Features

  * Change the temperature and pressure for log&#8201;*K* calculations.
  * Update parameters for activity and osmotic coefficients (Debye-Hückel, H<sub>2</sub>O, CO<sub>2</sub>) for specified *T* and *P*.
  * Add new species to the output file.
  * List data references in the output file (requires CHNOSZ version > 1.3.6).
  * Make plots comparing log&#8201;*K* values from two data files.
  * Update OBIGT using thermodynamic parameters (Δ*G*°, *S*° and <i>C<sub>p</sub></i>°) fitted to log&#8201;*K* values from a GWB file.
    * Enables addition of species from a GWB file to calculations and diagrams in CHNOSZ.

## Installation

```R
install.packages("remotes")
remotes::install_github("jedick/logKcalc", build_vignettes = TRUE)
```

## Documentation

See the [vignettes](http://chnosz.net/#logKcalc) for detailed usage examples.
