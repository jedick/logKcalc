# logKcalc

This R package calculates equilibrium constants for thermodynamic data files used in [The Geochemist's Workbench®](https://www.gwb.com).
It reads reactions from an input GWB data file, calculates equilibrium constants (log*K*) using the [OBIGT database in CHNOSZ](http://chnosz.net/vignettes/obigt.html), and writes them to an output GWB file.

## Installation

```R
install.packages("remotes")
remotes::install_github("jedick/logKcalc", build_vignettes = TRUE)
```

## Features

  * Ability to change the temperature and pressure for log*K* calculations.
  * List data references in the output file (requires CHNOSZ version > 1.3.6).
  * Make plots comparing log*K* values from two data files.

## TODO (planned features)

  * Update activity coefficient model blocks (Debye-Hückel, H<sub>2</sub>O, CO<sub>2</sub>) for specified *T* and *P*.
  * Add new species; currently, only species in the input file are included.
