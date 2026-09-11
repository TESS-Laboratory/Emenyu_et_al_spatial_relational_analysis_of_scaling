# Spatial and relational dynamics of agroforestry adoption: integrating spatial and social network approaches in Uganda

This repository contains the data and code used to reproduce the analyses, figures and tables presented in:

Emenyu, A.P., Cunliffe, A. M., Jasny, L. & Powell, T.W.R. (2026) *Spatial and relational dynamics of agroforestry adoption: integrating spatial and social network approaches in Uganda*

A permanent version of this repository is archived at [https://zenodo.org/badge/1306576381.svg)](https://doi.org/10.5281/zenodo.22686984)

 Use of this code is licensed under [![Licence: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](...)
 
Contact Author: Antony P. Emenyu Email: ae474@exeter.ac.uk / emanyoph@gmail.com
 
 ## Project overview

 This project explores the spatial relational dynamics associated with TIST tree planting patterns in Uganda. This repository provides the analytical tools and a cleaned dataset to reproduce the analysis.

 ### Repository contents
 This repository has two key folders. The data folder and Scripts folder.

 ## Data folder
This folder contains the cleaned-pseudonymised dataset used in the analysis.
## TISTDat folder
This folder contains two .csv files. the Cleaned_anonymised_.. is the primary dataset covering the entire country. The BushSoroti_Cleaned... dataset however is a subset of the later with the subcounty columns added for comparison analysis between the two subcounties.
## UGshapefiles and Uganda_Forest_Reserves
These two folders contain shape files  for Uganda administrative boundaries and forest reserves respectively.

## Scripts folder
Contains all the scripts used in the analysis- numbers in execution sequence from 1-5. When executing the scripts using the datasets in the Data folder, I recommend starting with script 2 onwards. Script 1 contains the data wrangling procedures and code used to generate the .csv filed in the data folder.

|Script | Description |
|*2. Descriptive_statistics* | National and cross-site comparative plots |
|*3. Landcover, planting...* | Landcover, planting suitability,forest proximity comparison |
|*4. Semi-variograms  .....* | Semi-variograms, Mixed effects models and associated visuals |
|*5. Multilayer network ...* | Exponential random graph models, planting proximity simulation |
