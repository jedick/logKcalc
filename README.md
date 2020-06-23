# logKcalc

This R package calculates equilibrium constants for thermodynamic data files used in [The Geochemist's Workbench®](https://www.gwb.com).
It reads reactions from an input GWB data file, calculates equilibrium constants (log*K*) for the reactions using the [OBIGT database in CHNOSZ](http://chnosz.net/vignettes/obigt.html), and writes them to an output GWB file.

## Installation

```R
install.packages("remotes")
remotes::install_github("jedick/logKcalc", build_vignettes = TRUE)
```

## Features

  * Change the temperature and pressure for log*K* calculations.
  * Update parameters for activity and osmotic coefficients (Debye-Hückel, H<sub>2</sub>O, CO<sub>2</sub>) for specified *T* and *P*.
  * Add new species to the output file.
  * List data references in the output file (requires CHNOSZ version > 1.3.6).
  * Make plots comparing log*K* values from two data files.
