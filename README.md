# Driving a protective single allelic variant of the mosquito FREP1 gene to combat malaria

This repository contains R code for modeling the dynamics of a protective single-allelic variant of the mosquito FREP1 gene across discrete generations.

# Repository Overview

This repository provides:

1. Scripts for simulating multi-generational non-overlapping dynamics of the FREP1 gene variant in mosquito populations under laboratory conditions.

2. Parameter estimation methods using Baesyan methods and experimental data.

3. Tools for visualizing simulation results.

# File Descriptions

<a href=https://github.com/lambsUSP/FREP1/blob/main/main_FREP1.R>main_FREP1.R</a>: Main script for running simulations and performing parameter estimations.

loading_data_experiment.R: Loads experimental data.

tMatrix_func.R: Defines the transition matrix for the model, capturing genotype dynamics.

maternal_deposition.R: Models the effects of maternal deposition on allele inheritance.

single_sim_model.R: Simulates a single generation step, capturing allele frequency changes.

experiment_sim_model.R: Simulates multiple non-overlapping generations.

obj_func.R: Computes the likelihood function for parameter estimation, comparing observed and simulated data.

metropolis_MCMC.R: Implements the Metropolis algorithm for parameter estimation using Markov Chain Monte Carlo (MCMC) methods.

plotting_data_experiment.R: Generates plots to visualize experimental data and simulation results.

# Getting Started

The main script main_FREP1.R loads the other functions and performe simulations and parameter estimations.

# Authors

Rodrigo M. Corder and John M. Marshall.
