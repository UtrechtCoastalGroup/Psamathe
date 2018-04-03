# Aeolus
Ground water - surface moisture - aeolian sand transport model for narrow sandy beaches

The Aeolus repository contains Matlab-scripts to predict meso-scale (seasonal to annual) aeolian input of beach sand into the foredune system. The Aeolus model has been motivated by two seminal papers:
(1) Bauer, B.O. and R.G.D. Davidson-Arnott, 2002. A general framework for modeling sediment supply to coastal dunes including wind angle, beach geometry, and fetch effects. Geomorphology, 49, 89-108.
(2) Delgado-Fernandez, I., 2011. Meso-scale modelling of aeolian sediment input to coastal dunes. Geomorphology, 2011, 230-243.

The Aeolus model has three modules:
(1) A beach groundwater model, based on the non-linear Boussinesq equation including the effect of wave runup.
(2) A surface-moisture model, which takes the output of the beach groundwater model and computes surface moisture content using the Van-Genuchten soil water retention curve.
(3) A fetch-based aeolian sand transport equation, in which the potential sand transport rate is computed according to the formulations of Hsu.  
The 'fetch' model includes the effect of spatially varying surface moisture on critical fetch and can handle the effect of the foredune on wind on the beach by a "beach wind angle" and a "foredune wind angle". Generally, the latter will be more alongshore than the former. The Hsu model has been 're-formulated' to include the dependence of aeolian sand transport on grainsize and the height of the wind sensor above the bed.

The groundwater - surface-moisture model has been described in (currently under review):
Brakenhoff, L.B., Y. Smit, J.J.A. Donker and G. Ruessink, under review. Tide-induced variability in beach surface moisture: observations and modelling. Earth Surface Processes and Landforms.

Funded by the Dutch Technology Foundation (STW) of the Netherlands Organisation for Scientific Research (NWO), Vici project #13709 awarded to Gerben Ruessink.

Gerben Ruessink - April 3, 2018
