# Aeolus: Model for computing wind-driven sand transport from the beach onto the foredune based on the fetch concept

### Description

The Aeolus repository contains Matlab-scripts to predict meso-scale (seasonal to annual) aeolian input of beach sand into the foredune system. The Aeolus model has been motivated by two seminal papers:
- Bauer, B.O. and R.G.D. Davidson-Arnott, 2002. A general framework for modeling sediment supply to coastal dunes including wind angle, beach geometry, and fetch effects. Geomorphology, 49, 89-108. https://doi.org/10.1016/S0169-555X(02)00165-4
- Delgado-Fernandez, I., 2011. Meso-scale modelling of aeolian sediment input to coastal dunes. Geomorphology, 2011, 230-243. https://doi.org/10.1016/j.geomorph.2011.04.001

The Aeolus model has three main modules:
- A beach groundwater module, based on the non-linear Boussinesq equation including the effect of wave runup;
- A surface moisture module that uses the Van Genunchten soil water retention curve to compute the spatio-temporal groundwater depth values from the first module into surface moisture;
- A fetch-based aeolian sand transport equation, in which the downwind increase in aeolian transport is computed based on the fetch concept; the critical fetch depends on the wind speed and the spatially varying surface moisture.

At this moment, only module 1 has been posted, directory Groundwater. Subdirectory example contains an example run file with corrosponding input data. 

The combined groundwater and surface-moisture modules have been described in:
- Brakenhoff, L.B., Y. Smit, J.J.A. Donker and G. Ruessink, 2019. Tide-induced variability in beach surface moisture: observations and modelling. Earth Surface Processes and Landforms, 44, 317-330. https://doi.org/10.1002/esp.4493 (Open Access)

A first version of the fetch model has been described in:
(1) Smit, Y., 2019. Surface moisture dynamics on a narrow coastal beach. PhD thesis, Utrecht University. ISSN 2211-4335.
(2) Hage, P., G. Ruessink, Z. van Aartrijk and J. Donker (under review). Using video monitoring to test a fetch-based aeolian sand transport model. Journal of Marine Science and Engineering.

We are currently working on a test of the model against a high-resolution, multi-annual data-set of foredune growth described in Ruessink et al. (2019), Data, 4, 73; doi:10.3390/data4020073. As soon as this paper has been submitted (probably mid-2020), the remaining code, with examples, will be made fully available.

Funded by the Dutch Technology Foundation (STW) of the Netherlands Organisation for Scientific Research (NWO), Vici project #13709 awarded to Gerben Ruessink.

Gerben Ruessink - February 4, 2020
