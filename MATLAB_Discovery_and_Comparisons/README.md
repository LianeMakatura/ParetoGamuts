# Pareto Gamuts: Exploring Optimal Designs Across Varying Contexts

Liane Makatura, Minghao Guo, Adriana Schulz, Justin Solomon, and Wojciech Matusik

SIGGRAPH 2021

This repository contains the code to accompany our SIGGRAPH 2021 paper, 
[Pareto Gamuts: Exploring Optimal Designs Across Varying Contexts](http://paretogamuts.csail.mit.edu). 
For more details, please refer to the paper.

## Overview

This folder contains (1) our discovery algorithm, (2) problem setups for all examples scenarios included in the paper, and (3) code for running comparisons against two state-of-the-art fixed context approaches (termed "NSGA" and "Schulz" -- see paper for details.)


## Pareto Gamut Discovery code

All Pareto gamut code is in the folder `ParetoGamuts_discovery`.
You must change the working directory in `AppDiscover.m`, line 12. Be sure to include the final `/` in the path.

While developing the code, we used the term "application variable" to describe what is known in the paper as the "context variable"; they should be treated interchangably here.


### Running the examples

All examples from Section 6.1 (ZDT1-6, Fourier benchmark) can be run from `examplesAnalytic.m`

All examples from Section 6.2 (turbine, bicycle, solar roof/gridshell, lamp, bicopter) can be run from `examplesMeshBased.m`


## Comparisons with NSGA

All code for NSGA is in the folder `NSGA_comparisons`. 

### Running the examples

Turbine and Lamp:  run `python comparison_stats.py`

To switch between the two problems, you will need to make some small changes within the file:
- line 53: change the nadir point for HV to 2 (Turbine) or 3(Lamp) dimensions
- line 165: change the problem class to Turbine or Lamp

Bikerocker: run `python runBike.py`

Bicopter: run `python runBicopter.py`

Solar Roof: run `python runGridShell.py`


## Comparisons with Schulz

All code for Schulz is in the folder `ParetoGamuts_discovery`.
You must change the working directory in `MappingFunction.m`, line 51. Be sure to include the final `/` in the path.

### Running examples

All examples can be run from `compareGamutSingleContext.m`. You will need to specify `scFunc_list` according to the descriptions in `MappingFunction.m`.

 