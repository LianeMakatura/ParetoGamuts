# Pareto Gamuts: Exploring Optimal Designs Across Varying Contexts
## SIGGRAPH 2021 Submission ID 153


This repository contains 
* our discovery algorithm, 

* problem setups for all examples scenarios included in the paper, and 

* code for running comparisons against two state-of-the-art fixed context approaches (termed "NSGA" and "Schulz" -- see paper for details.)


# Pareto Gamut Discovery code

## Before getting started
All Pareto gamut code is in the folder `cleanUp_2`.
You must change the working directory in `AppDiscover.m`, line 12. Be sure to include the final `/` in the path.


## Analytic Examples from Section 6.1 (ZDT1-6, Fourier benchmark)
All examples can be run from `examplesAnalytic.m`

## Engineering examples from Section 6.2 (turbine, bicycle, solar roof/gridshell, lamp, bicopter)
All examples can be run from `examplesMeshBased.m`

## Background Information
While developing the code, we used the term "application variable" to describe what is known in the paper as the "context variable"; they should be treated interchangably here.

# Comparisons with NSGA
## Before getting started
All code for NSGA is in the folder `NSGA_comparisons`. 
## Running the examples
Turbine and Lamp:  run `python comparison_stats.py`

Need some little changes of the file in:
- line 53: change the nadir point for HV to 2 (Turbine) or 3(Lamp) dimensions
- line 165: change the problem class to Turbine or Lamp

Bikerocker: run `python runBike.py`

Bicopter: run `python runBicopter.py`

Solar Roof: run `python runGridShell.py`


# Comparisons with Schulz

## Before getting started 
All code for Schulz is in the folder `cleanUp_2`.
You must change the working directory in `MappingFunction.m`, line 51. Be sure to include the final `/` in the path.

## Running examples
All examples can be run from `compareGamutSingleContext.m`. You will need to specify `scFunc_list` according to the descriptions in `MappingFunction.m`.

