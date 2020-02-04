# Aeolus
Ground water - surface moisture - aeolian sand transport model for narrow sandy beaches

The Aeolus repository contains Matlab-scripts to predict meso-scale (seasonal to annual) aeolian input of beach sand into the foredune system. The Aeolus model has been motivated by two seminal papers:
(1) Bauer, B.O. and R.G.D. Davidson-Arnott, 2002. A general framework for modeling sediment supply to coastal dunes including wind angle, beach geometry, and fetch effects. Geomorphology, 49, 89-108.
(2) Delgado-Fernandez, I., 2011. Meso-scale modelling of aeolian sediment input to coastal dunes. Geomorphology, 2011, 230-243.

The Aeolus model has two modules:
(1) A beach groundwater - surface moisture model, based on (groundwater) the non-linear Boussinesq equation including the effect of wave runup and (surface moisture) a Van Genuchten soil water retention curve;
(2) A fetch-based aeolian sand transport equation, in which the potential sand transport rate is computed according to the formulations of Hsu, De Kok or Lettau-Lettau.  

The 'fetch' model includes the effect of spatially varying surface moisture on critical fetch and can handle the effect of the foredune on wind on the beach by a "beach wind angle" and a "foredune wind angle". Generally, the latter will be more alongshore than the former (De Winter et al., 2020).

At this moment, only module 1 has been posted, directory Groundwater. Subdirectory example contains an example run file with corrosponding input data. 

The groundwater - surface-moisture model has been described in:
Brakenhoff, L.B., Y. Smit, J.J.A. Donker and G. Ruessink, 2019. Tide-induced variability in beach surface moisture: observations and modelling. Earth Surface Processes and Landforms, 44, 317-330.

A first version of the fetch model has been described in:
Smit, Y., 2019. Surface moisture dynamics on a narrow coastal beach. PhD thesis, Utrecht University. ISSN 2211-4335.
Hage, P., G. Ruessink, Z. van Aartrijk and J. Donker (under review). Using video monitoring to test a fetch-based aeolian sand transport model. Journal of Marine Science and Engineering.

We are currently working on a test of the model against a high-resolution, multi-annual data-set of foredune growth described in Ruessink et al. (2019), Data, 4, 73; doi:10.3390/data4020073. As soon as this paper has been submitted (probably mid-2020), the remaining code, with examples, will be made open-access.

Funded by the Dutch Technology Foundation (STW) of the Netherlands Organisation for Scientific Research (NWO), Vici project #13709 awarded to Gerben Ruessink.

Gerben Ruessink - February 4, 2020
