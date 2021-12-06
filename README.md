# Psamathe: Modelling wind-driven sand transport from the beach using the fetch concept

### Description

The Psamathe repository contains Matlab-scripts to predict meso-scale (seasonal to annual) aeolian sand supply from the beach to the foredune. The Psamathe model, named after the Greek goddess of sand beaches, was motivated by two seminal papers:
- Bauer, B.O. and R.G.D. Davidson-Arnott, 2003. A general framework for modeling sediment supply to coastal dunes including wind angle, beach geometry, and fetch effects. Geomorphology, 49, 89-108. https://doi.org/10.1016/S0169-555X(02)00165-4
- Delgado-Fernandez, I., 2011. Meso-scale modelling of aeolian sediment input to coastal dunes. Geomorphology, 2011, 230-243. https://doi.org/10.1016/j.geomorph.2011.04.001

The Psamathe model has three main modules:
- A beach groundwater module, based on the non-linear Boussinesq equation including the effect of wave runup;
- A surface moisture module that uses the Van Genunchten soil water retention curve to compute the spatio-temporal groundwater depth values from the first module into surface moisture;
- A fetch-based aeolian sand transport equation, in which the downwind increase in aeolian transport is computed based on the fetch concept; the critical fetch depends on the wind speed and the spatially varying surface moisture.

### Groundwater and surface moisture

Modules 1 and 2 are posted in the directory Groundwater. Its subdirectory "example" contains an example run file with corrosponding input data. 

The combined groundwater and surface-moisture modules have been described in:
- Brakenhoff, L.B., Y. Smit, J.J.A. Donker and G. Ruessink, 2019. Tide-induced variability in beach surface moisture: observations and modelling. Earth Surface Processes and Landforms, 44, 317-330. https://doi.org/10.1002/esp.4493 (Open Access)

### Fetch

The third module is posted in the directory Fetch. Its subdirectory "example" contains an example run file with corresponding input data.

The fetch module has seen various versions, as witnessed in the following papers: 
- Hage, P., G. Ruessink, Z. van Aartrijk and J. Donker, 2020. Using video monitoring to test a fetch-based aeolian sand transport model. Journal of Marine Science and Engineering. https://doi.org/10.3390/jmse8020110 (Open Access)
- Tuijnman, J., J.J.A. Donker, C.S. Schwarz and G. Ruessink, 2020. Consequences of a storm surge for aeolian sand transport on a low-gradient beach. Journal of Marine Science and Engineering. https://doi.org/10.3390/jmse8080584 (Open Access)

In both papers the model was called Aeolus. But there appeared to be other Aeolus models, so Aeolus became Psamathe.

The full model, including the most recent modifications to the fetch module, can be found in:

- Ruessink, G., G. Sterk, Y. Smit, W. de Winter, P. Hage, J.J.A. Donker and B. Arens, under review. Predicting monthly to multi-annual foredune growth at a narrow beach. Earth Surface Processes and Landforms. (Will be updated once accepted for publication)

Funded by the Dutch Technology Foundation (STW) of the Netherlands Organisation for Scientific Research (NWO), Vici project #13709.

Gerben Ruessink - December 6, 2021
