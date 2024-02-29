# Spatstat_Yiming
## Introduction
This folder contains datasets and R-code corresponding to the report "Estimating animal resource selection from telemetry data using space-time point process models" written by Yiming Zhao. A space-time model is constructed and applied to the datasets for resource selection function (RSF) analysis of a brown bear in southeast Alaska.

The folder is structured as follows:
- File bear_track.csv contains the telemetry locations of the bear.
- File bear_habitat.csv contains two environmental covariates within the bear habitat.
- File gray.colors.rev.R provides a function to generate grayscale color vector when plotting.
- File draw.contour.R provides a function to draw contour lines for datasets.
- File draw.contour.2.R provides a function similar to the function in draw.contour.R.
- File stpp_rsf_helper.R provides functions for the construction of the quadrature points and the space-time point process model.
- File run_bear_analysis.R analyzes RSF of the bear data with the space-time model.
- File plot.R plots the bear data and the model fitting results.
- File setup_session.R installs all needed R-packages.

The code was tested on R version 4.2.1 on windows 10.

## Setup
First, please install all necessary packages in file setup_session.R.

Then, please download the datasets:
- bear_track.csv
- bear_habitat.csv

Afterwards, please download and run the following files: 
- gray.colors.rev.R
- draw.contour.R
- draw.contour.2.R

In order to reproduce the results, please download the rest R files and run them in the following order:
- stpp_rsf_helper.R
- run_bear_analysis.R (estimated runtime: 20min)
- plots.R
The file stpp_rsf_helper.R provides necessary functions. It is then sourced by run_bear_analysis.R to apply the models to the datasets provided and all the results reproduced are stored in a new file bear_stpp_results.RData. The file plots.R sources the results in bear_stpp_results.RData for explanation of the datasets and visualization the model fitting results.

## References to the datasets
- All the datasets and code: Johnson, D.S., corresponding to Chapter 4 of the book: Animal movement: statistical models for telemetry data (2017), CRC press. secs. 4.7.3-4.7.4. Everything was accessed via email on December 13th, 2023.

