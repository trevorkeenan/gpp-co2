# GPP-CO2
## Overview
This repository contains analysis and plotting scripts to reproduce the emergent constraint analysis presented in
Keenan et al. 2021: A constraint on historic growth in global photosynthesis due to increasing CO2.
 
Nature - https://www.nature.com/articles/s41586-021-04096-9

Full information on the methods used in this study are attached to this paper and are available
online at the link provided above; this includes information about datasets used, as well as the motivation and reasoning
behind the analysis.

The scripts, written in Matlab, are:
## A. varianceNormalization.m
This script derives the relationship between Sland and Beta^{GPP}, and performs variance normalization 
to extract the partial response. 

## B. calc_EC_andPlot.m
This code, called by 'varianceNormalization.m' uses the emergent constraint between 
the partial response of Sland to Beta^{GPP} across models to derive the constrained Beta^{GPP} 

## Running the code
Running A_varianceNormalization.m will produce the following figures reported in Keenan et al. 2021: \
*Fig. 1a-d  \
*ED Figure 1  \
*ED Figure 2  \
*ED Figure 3  \
*ED Figure 6  \
 
A_varianceNormalization.m calls B_calc_EC_andPlot.m

## Folder structure and contents
 ./figures/emergent contains the figures produced 
if the save_figures flag is set to 1.

./TRENDYv6_derived contains derived output from the TRENDY model simulations. 
TRENDY model simulations are not publically available but can be obtained through request to Prof. Sitch (S.A.Sitch@exeter.ac.uk)

./dataIntermediates contains output from the scripts included here, 
extracted from the data contained in ./TRENDYv6_derived

./functions contains plotting code and the prediction error code.


