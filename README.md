# Reduced-Dimensional Gravity Inversion

MATLAB implementation of the inversion framework presented in

Search Space Equivalence and Reduced-Dimensional Parametrizations for 3D Gravity Inversion Using RR-GPSO

developed by the GPINV Research Group, University of Oviedo.



# Overview

This repository contains the complete MATLAB implementation used to perform the numerical experiments described in the accompanying manuscript.

The software implements a stochastic 3D gravity inversion framework based on Real-Representation Particle Swarm Optimization (RR-GPSO) and compares different reduced-dimensional geometric parameterizations under physically comparable search-space constraints.

The package includes:

* synthetic inversion experiments;
* real-data inversion experiments;
* forward modeling routines;
* RR-GPSO implementation;
* statistical analysis scripts;
* real gravity datasets.

The objective is to evaluate how different geometric parameterizations influence inversion accuracy, convergence behavior, model complexity, and computational cost while maintaining physically comparable admissible search domains.



# Repository Structure

CaliforniaISO.xlsx                Real gravity dataset (Clear Lake)
Volcan-1.xlsx                     Real gravity dataset (Nirano)

MAIN\_PSO\_SINTETICO.m              Synthetic experiments
MAIN\_PSO\_SINTETICO\_REC.m          Synthetic experiments (rectangular prisms)
MAIN\_REAL\_RECT.m                  Real-data inversion

pso\_grav3D.m                      RR-GPSO algorithm
PSO\_options.m                     Optimization parameters
initialpop.m                      Initial population generation

GravPrismaCuadrado.m              Forward model (square prism)
GravPrismaRectangular.m           Forward model (rectangular prism)
GravElipProlato.m                 Forward model (prolate ellipsoid)
GravElipsoideProlatoY.m           Auxiliary ellipsoid routine
GravNubePuntos.m                  Point-cloud forward model

DiscretizaPrisma.m
DiscretizaPrismaRectangular.m
DiscretizaElipsoideProlato.m      Geometry discretization routines

fcost\_\*.m                         Objective functions

problema\_directo\_prolato.m        Synthetic forward modeling

estadisticas\_\*.m                  Statistical analysis scripts

dibuja\_anomali.m                  Plotting utilities

trata\_modeli.m
trata\_modelip.m                   Post-processing routines

posterior.m                       Posterior model analysis

cloud\_RR.mat                      Auxiliary data



# Software Requirements

* MATLAB
* The code is self-contained and does not depend on external libraries.



# Included Datasets

Two real gravity datasets are included.

## CaliforniaISO.xlsx

Clear Lake gravity anomaly used as a benchmark real-data example.

## Volcan-1.xlsx

Microgravity observations acquired over the Nirano mud volcano field.

Both datasets are directly read by the inversion scripts.



# Running the Code

## Synthetic experiments

Execute

matlab
MAIN\_PSO\_SINTETICO
or
matlab
MAIN\_PSO\_SINTETICO\_REC
depending on the selected parameterization.

These scripts generate synthetic gravity data, perform RR-GPSO inversion, and save the inversion results.



## Real-data inversion

Run
matlab
MAIN\_REAL\_RECT

The script reads the selected dataset, performs the inversion, and stores the results.



## Statistical analysis

After completing multiple independent inversion runs, execute the corresponding statistical scripts:
estadisticas\_sintetico.m
estadisticas\_prismas\_Cali.m
estadisticas\_prolato\_Cali.m
Estadisticas\_Nirano.m

These scripts compute descriptive statistics of the final relative error and generate the numerical summaries reported in the manuscript.



# Methodology

The inversion framework is based on

* Real-Representation Particle Swarm Optimization (RR-GPSO)
* physically comparable search-space constraints
* reduced-dimensional geometric parameterizations
* gravity forward modeling
* stochastic optimization using multiple independent runs.

The framework allows comparison of different geometric representations while maintaining comparable physical constraints.



# Expected Outputs

The inversion scripts generate

* recovered models;
* predicted gravity anomalies;
* convergence histories;
* execution times;
* MATLAB result files (.mat).

The statistical scripts compute

* minimum error;
* median error;
* mean error;
* standard deviation;
* interquartile range;
* Wilcoxon rank-sum tests.



# Reproducibility

The repository contains all files required to reproduce the experiments described in the accompanying manuscript.

No external datasets or third-party software are required.

Because RR-GPSO is a stochastic optimization algorithm, individual runs may differ slightly. The statistical analyses reported in the paper are based on multiple independent realizations.



# Citation

If you use this software in academic work, please cite the accompanying publication.

Search Space Equivalence and Reduced-Dimensional Parametrizations for 3D Gravity Inversion Using RR-GPSO



# Authors

GPINV Research Group

Department of Mathematics

University of Oviedo

Spain



# License

This project is distributed under the MIT License.

Copyright (c) GPINV Research Group

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, subject to inclusion of the above copyright notice and this permission notice.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED.

