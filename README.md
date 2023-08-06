# Biophysical model of a cell-free gluconate responsive biosensor

### Introduction
This repository contains code written in the [Julia](https://www.julialang.org) programming language for the cell-free gluconate responsive biosensor described in the publication: 

[Adhikari et al. (2023) Modeling and Analysis of a Cell-Free Gluconate Responsive Biosensor (bioRxiv)](https://www.biorxiv.org/content/10.1101/2023.01.10.523462v1.full)


### Installation and Requirements

To get the model codes, you can download this model repository as a zip file, clone or pull it using the command:

    git pull https://github.com/varnerlab/GntR-model-code-publication.git

or

    git clone https://github.com/varnerlab/GntR-model-code-publication.git

The `src` directory contains the code for the model, the `data` directory contains the experimental data, the `simulation` directory contains the simulated model files, `config` contains the configuration files for the model, and the `figs` directory contains the figures generated by the model.

To install the required Julia packages, run the following command in the Julia REPL:

    julia> using Pkg
    julia> Pkg.activate(".")
    julia> Pkg.instantiate()

### Scripts
Script | Description
---: | ---
`RUN_ALL.jl` | Runs all the required scripts for the model
`parameter_estimation.jl` | Solves the model equations for the ensemble of parameter sets for the test case, 10mM gluconate. Saves solutions in the `simulations/poets_ensemble` directory
`gluconate_dynamics.jl` | Uses the parameter ensemble generated by the `parameter_estimation.jl` script to simulate the TX/TL dynamics for varying gluconate concentrations. Saves the solutions in the ``simulated/gluconate_dynamics`` directory. Subdirectories within `gluconate_dynamics` hold solutions for each specified gluconate concentration
`dose_response_ensemble.jl` | Solves the model equations to generate the dose-response curve. Results are saved in the `simulated/dose_response_simulations` directory
`sensitivity_run.jl` | Conducts [Morris sensitivity analysis](https://doi.org/10.2307%2F1269043) on the model and saves the results in the `simulated/sensitivity_results` directory
`plot_*` | Plot the TX/TL dynamics from the simulations generated by the `gluconate_dynamics.jl` script, the dose-response curve from the simulations
