# FLORENCE

FLORENCE is an R package that provides functions to compute expected values of the site frequency spectrum of somatic variants generated by dividing tissues that are either expanding or homeostatic. Modalities for both, neutrally evolving tissues and tissues with selection of a single clone are provided. The computations assume deterministic cell numbers in the tissue overall but account for the stochasticity of clonal growth at the level of individual clones, using clone size distributions generated by linear birth-death processes. 

FLORENCE also provides modalities to simulate these processes stochastically, using gillespie's algoirthm. In addition, it allows simulation of VAF distributions generated by bulk whole genome sequencing or by pseudo-bulk data from single-cell WGS. FLORENCE was developed as part of the study "Detecting and quantifying clonal selection in somatic mosaics". Please refer to the github repository [clonal_hemopoiesis](https://github.com/VerenaK90/clonal_hemopoiesis) for further information.

## Documentation

Function documentation is available in the accompanying manual (FLORENCE_0.0.0.9000.pdf)[FLORENCE_0.0.0.9000.pdf].

## Demo 

An exemplary guide on how to use FLORENCE to model variant accumulation and to do parameter estimation in conjunction with pyABC is given in the accompanying vignette. Re-running this example cases should take 10-15 minutes.

## System requirements

### Hardware requirements

FLORENCE runs on a standard computer.

### Software requirements

#### OS Requirements

This package is supported for macOS and Linux. The package has been tested on the following systems:

macOS: Monterey (12.6.3)
Linux: CentOS Linux 7 (Core)

#### R dependencies

FLORENCE has been tested on R v4.2.0 and v4.2.1 and requires installation of the packages deSolve (v1.34), ape (5.6.2), phangorn (v2.10.0), phytools (v1.2.0).

FLORENCE can be run in conjunction with python3 and pyABC (≥0.12.6) to estimate parameters from whole genome sequencing data.

## Version

0.0.0.9000

## Citation

Körber et al., Detecting and quantifying clonal selection in somatic mosaics.

## Installation

library(devtools)
devtools::install.github("VerenaK90/FLORENCE")

If you want to install the vignette, run

devtools::install.github("VerenaK90/FLORENCE", build_vignettes=TRUE). The vignette can then be viewed via vignette("Vignette", "FLORENCE"). The package including its vignette can be installed within 10-15 minutes given that python, R and pyABC have been installed.

## License

FLORENCE is run under an MIT License.

## Contact

Verena Körber (v.koerber@dkfz.de)
