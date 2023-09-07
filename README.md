# SCI&Phi;ER

SCI&Phi;ER is an R package that provides functions to compute expected values of the site frequency spectrum of somatic variants generated by dividing tissues that are either expanding or homeostatic. Modalities for both, neutrally evolving tissues and tissues with selection of a single clone are provided. The computations assume deterministic cell numbers in the tissue overall but account for the stochasticity of clonal growth at the level of individual clones, using clone size distributions generated by linear birth-death processes. 

SCI&Phi;ER also provides modalities to simulate these processes stochastically, using gillespie's algoirthm. In addition, it allows simulation of VAF distributions generated by bulk whole genome sequencing or by pseudo-bulk data from single-cell WGS. SCI&Phi;ER was developed as part of the study "Detecting and quantifying clonal selection in somatic mosaics" by Verena Körber et al. Please refer to the github repository [Clonal_hematopoiesis](https://github.com/VerenaK90/clonal_hematopoiesis) for further information.

## Documentation

Function documentation is available in the accompanying manual [SCIFER_0.0.0.9000.pdf](SCIFER_0.0.0.9000.pdf).

## Demo 

An exemplary guide on how to use SCI&Phi;ER to model variant accumulation and to do parameter estimation in conjunction with pyABC is given in the accompanying vignette. Re-running this example cases should take 10-15 minutes.

## System requirements

### Hardware requirements

SCI&Phi;ER runs on a standard computer.

### Software requirements

#### OS Requirements

This package is supported for macOS and Linux. The package has been tested on the following systems:

macOS: Monterey (12.6.3)
Linux: CentOS Linux 7 (Core)

#### R dependencies

SCI&Phi;ER has been tested on R v4.2.0 and v4.2.1 and requires installation of the packages deSolve (v1.34), ape (5.6.2), phangorn (v2.10.0), phytools (v1.2.0).

SCI&Phi;ER can be run in conjunction with python3 and pyABC (?0.12.6) to estimate parameters from whole genome sequencing data.

## Version

0.0.0.9000

## Citation

Körber et al., Detecting and quantifying clonal selection in somatic mosaics.

## Installation

library(devtools)
devtools::install_github("VerenaK90/SCIFER/tree/paper", ref="paper")

If you want to install the vignette, run

devtools::install.github("VerenaK90/SCIFER/tree/paper", ref="paper", build_vignettes=TRUE). The vignette can then be viewed via vignette("Vignette", "SCIFER"). The package including its vignette can be installed within 10-15 minutes given that python, R and pyABC have been installed.

## License

SCI&Phi;ER is run under an [MIT License](https://web.archive.org/web/20160411224647/https://opensource.org/licenses/MIT).

## Contact

Verena Körber (v.koerber@dkfz.de)
