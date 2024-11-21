# Driving a protective single allelic variant of the mosquito FREP1 gene to combat malaria

This repository contains R code for modeling the dynamics of a protective single-allelic variant of the mosquito FREP1 gene across discrete generations.

# Repository Overview

This repository provides:

1. Scripts for simulating multi-generational non-overlapping dynamics of the FREP1 gene variant in mosquito populations under laboratory conditions.

2. Parameter estimation methods using Baesyan methods and experimental data.

3. Tools for visualizing simulation results.

# File Descriptions

<a href=https://github.com/lambsUSP/FREP1/blob/main/main_FREP1.R>main_FREP1.R</a>: Main script for running simulations and performing parameter estimations.

<a href=https://github.com/lambsUSP/FREP1/blob/main/loading_data_experiment.R>loading_data_experiment.R</a>: Loads experimental data.

<a href=https://github.com/lambsUSP/FREP1/blob/main/tMatrix_func.R>tMatrix_func.R</a>: Defines the transition matrix for the model, capturing genotype dynamics.

<a href=https://github.com/lambsUSP/FREP1/blob/main/maternal_deposition.R>maternal_deposition.R</a>: Models the effects of maternal deposition on allele inheritance.

<a href=https://github.com/lambsUSP/FREP1/blob/main/single_sim_model.R>single_sim_model.R</a>: Simulates a single generation step, capturing allele frequency changes.

<a href=https://github.com/lambsUSP/FREP1/blob/main/experiment_sim_model.R>experiment_sim_model.R</a>: Simulates multiple non-overlapping generations.

<a href=https://github.com/lambsUSP/FREP1/blob/main/mobj_func.R>obj_func.R</a>: Computes the likelihood function for parameter estimation, comparing observed and simulated data.

<a href=https://github.com/lambsUSP/FREP1/blob/main/metropolis_MCMC.R>metropolis_MCMC.R</a>: Implements the Metropolis algorithm for parameter estimation using Markov Chain Monte Carlo (MCMC) methods.

<a href=https://github.com/lambsUSP/FREP1/blob/main/plotting_data_experiment.R>plotting_data_experiment.R</a>: Generates plots to visualize experimental data and simulation results.

# Getting Started

The main script <a href=https://github.com/lambsUSP/FREP1/blob/main/main_FREP1.R>main_FREP1.R</a> loads the other functions and performe simulations and parameter estimations.

# Authors

<a href=https://orcid.org/0000-0002-8929-2200>Rodrigo M. Corder</a> and <a href=https://orcid.org/0000-0003-0603-7341>John M. Marshall</a>.
