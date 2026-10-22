# Geometric and demographic effects explain contrasting fragmentation-biodiversity relationships across scales

<!-- badges: start -->
<!-- badges: end -->

*This repository contains the code and data for our study on how habitat fragmentation affects biodiversity across spatial scales. Using a process-based metacommunity model, we explore how geometric and demographic effects shape fragmentation–biodiversity relationships. Below, you’ll find details on how to reproduce our results and run the model.*

## Manuscript

Gelber*, S., Blowes, S. A., Chase, J. M., Huth, A., Schurr, F. M., Tietjen, B., Zeller, J. W. and May, F. 2025. Geometric and demographic effects explain contrasting fragmentation-biodiversity relationships across scales. – Oikos, DOI: 10.1002/oik.10778.

*stav.gelber@fu-berlin.de

- [Read the preprint](https://www.biorxiv.org/content/10.1101/2024.02.01.577731v1)
- [![DOI](https://zenodo.org/badge/731115734.svg)](https://doi.org/10.5281/zenodo.14944795)

## Abstract

There is consensus that habitat loss is a major driver of biodiversity loss, while the effects of fragmentation, given a constant total habitat amount, are still debated. Here, we use a process-based metacommunity model to show how scale- and context-dependent fragmentation–-biodiversity relationships can emerge from the interplay of two types of fragmentation effects – geometric and demographic. Geometric effects arise from the spatial distributions of species and landscape modification, whereas demographic effects reflect long-term changes in species demographic rates following landscape modification. Our spatial model considers sessile individuals in a heterogeneous landscape and dynamically simulates the processes of species reproduction, dispersal, competition, mortality, and immigration. We introduce a novel approach to partition geometric and demographic fragmentation effects that is based on model outputs directly after landscape modification and after a phase of community dynamics in the modified landscape. In detailed simulation experiments, we assessed how key ecological processes and factors, such as dispersal, habitat heterogeneity, and edge effects, influence geometric, demographic, and net fragmentation effects across spatial scales. We found that increasing intraspecific aggregation due to short dispersal and/or environmental autocorrelation increased positive geometric fragmentation effects at the landscape scale. In our model, negative demographic fragmentation effects emerged at the local and landscape scale due to high dispersal mortality in the matrix and due to negative edge effects. We showed that the model can simultaneously predict positive fragmentation-biodiversity relationships at the local scale and negative relationships at the landscape scale as well as context-dependent variation of these relationships at the landscape scale. We conclude that the framework of geometric and demographic effects can reconcile previous apparently conflicting results and hopefully unlock and advance the debate on biodiversity changes in modified landscapes.

## Structure of the repository

-   `data-raw/` – Contains raw model simulation data used to generate figures. Each dataset includes corresponding model parameters.
-   `Model/` – Contains source code for running simulations, along with R scripts for executing the model locally or on a computing cluster.
-   `R/` – R scripts for generating figures from raw data.
  
## How to reproduce the figures

1. Clone this repository to your computer.

2. Unzip the output files under `data-raw/model_output/fig_s3a` and `data-raw/model_output/fig_s3b`

3. Run the script `R/generate_figures.R`

## How to run the model

1. Adjust the parameter file (`Model/parameters.R`). To reproduce previous runs, you can adjust it according to the static parameter file in each of the data subfolders under `data-raw/model_output`.

2. Execute the model using `Model/cluster_model_run.R` for cluster runs or `Model/run_model.R` for local runs.

3. Output files will be written to a new folder (`Model/outputs`) and will include a general and sample-scale output file for each model repetition, as well as separate files for the model's static and varying parameters.
