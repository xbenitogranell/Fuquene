# Fuquene
Multiproxy analysis of the 250-yr Lake Fuquene record (submitted Communications Earth and Environment).
This repository contains the necessary data and code to reproduce the statistical analyses and generate the main plots. 

## Study area and datasets description
+ **Lake Fuquene**: Pleistocene in origin, morrained-dammed tropical lake located in the Eastern Cordiellera of Colombian Andes. Water depth=2-6 m

## Two sediment cores, composite of 250 kyr diatom, geochemical and pollen record
+ **F9**: 60-m long, sampled for diatom, pollen and geochemistry
+ **F7**: sampled for diatoms and pollen (littoral)

## Main numerical steps:
+ Clean data and estimate Principal Curves for each assemblage (pollen and diatoms)
+ Plot stratigraphical floristic changes over time to identify periods of change
+ Interpolate PrC from coarser, explanatory variables (aquatic and terrestrial PrC) to lower resolution, responsive variables (diatom abundances) resolution 
+ Functional diversity analyses: how changes in lake habitat size did influence niche occupation by functionally similar diatom assemblages?
+ Distributed Hierarchical Generalized Additive Models (DHGAM)

## Folders
+ scripts: series of analyses arranged in a sequential order, including functions
+ data: individual files to read