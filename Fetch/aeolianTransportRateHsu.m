function qHsu = aeolianTransportRateHsu(U,par)
%AEOLIANTRANSPORTRATEHSU Hsu equation
% This function computes the potential aeolian sand transport rate given wind
% velocity, grain size and a factor that related the mean wind velocity at
% a height to the atmospheric shear velocity following the approach of Hsu (1974),
% and his earlier work.
%
% Additional comments:
% Hsu was rather confusing, in my humble opinion, on the units applied and
% I believe that in the definition of the Froude number Hsu was wrong.
% There D is said to be in mm, while U* and g have cm. So, D should be in
% cm too, to have the Froude number in a consistent, non-dimensional
% manner.
%
% The Hsu equation does not contain a threshold of motion.
% Delgado-Fernandez, who used the Hsu equation in her meso-scale modelling
% study, "filtered" the input wind data by setting a threshold wind
% velocity below which q was set to 0. This approach can be chosen here as
% well (see INPUT below), although not with a wind velocity but with an
% internally computed saltationFluidThreshold.
%
% INPUT
%  U, wind velocity, N values (row or column) of wind velocity (m/s)
%  par, a structure with parameter values. Of relevance here are the
%      fields g, a, D50 and thresholdWind. This threshold is a logical value. If
%      1 (i.e., use a threshold), then AN, rhoS, rhoA, D50 and gamma are
%      additionally required. 
%
% OUTPUT
%  qHsu, aeolian transport rate for each U, kg/m/s
%
% v1.0, Gerben Ruessink, November 17, 2017
% v1.1, Gerben Ruessink, April 2, 2019. Modifications: (1) code rewritten to be
%    more consistent with Aeolus model and other potential transport
%    equations (use of par); (2) implemented thresholdWind. This implies that, in
%    contrast, to earlier versions of the fetch model, there is no need
%    anymore to specify a minimum wind speed (as in Delgado-Fernandez's
%    work, for example).
% v1.2, Gerben Ruessink, May 13, 2019. Included rain effect (in uStar_it).
%
% REFERENCES
% Hsu, S.A., 1971. Wind stress criteria in eolian sand transport. Journal
%   of Geophysical Research, 76, 8684-8686.
% Hsu, S.A., 1974. Computing eolian sand transport from routine weather
%   data. Proc. 14th Int. Conf. on Coastal Engineering, 1619-1626.
% Delgado-Fernandez, I., 2011. Meso-scale modelling of aeolian sediment
%   input to coastal dunes. Geomorphology, 130, 230-243.

% error checking
fieldPresence = isfield(par,{'a','D50','g','thresholdWind','CRain','rainPotential'});
if sum(fieldPresence) ~= 6
    error('One or more parameters missing in aeolianTransportRateHsu.m. Required: a, g, D50, thresholdWind, rainPotential and CRain.');
end

% saltation fluid trenshold based on the method of Shao and Lu (2000), if
% needed
if par.thresholdWind
    uStar_it = 100*saltationFluidThreshold(par);  % in cm/s for consistency below
    
    % include effect of rain if so desired based on Arens (1996)
    if par.rainPotential
       uStar_it = uStar_it*(1+par.CRain);
    end
    
else
    uStar_it = 0;                                 % no threshold
end

% The equation requires diameter in mm (Dmm) and in cm (Dcm). The input
% here is in m. Thus:
Dmm = par.D50*1000;
Dcm = par.D50*100;

% qHsu is in g/cm/s, we would like to have qHsu in kg/m/s. See also
% Davidson-Arnott and Law (1996).
unitConversion = 0.1;

% Here is the Eolian transport coefficient K, in g/cm/s, with D = Dmm. The
% 10^-4 follows from the 10^4 along the y-axis of Figure 3, Hsu (1971).
K = exp(-0.47+4.97*Dmm)*10^-4;

% Froude number
g = 100*par.g;                  % cm/s2 as in Hsu papers
uStar = 100*par.a*U;            % uStar in cm/s, U in m/s as in Hsu papers
Fr = uStar / sqrt(g*Dcm); 

% Transport in kg/m/s
qHsu = unitConversion*K*Fr.^3;
qHsu(uStar < uStar_it) = 0;

% This results in: qHsu = 1.14x10^-5*U^3 (kg/m/s) and qHsu = 1.14x10^-4*U^3
% (g/cm/s). Correct. It also implies I interpreted Dmm and Dcm correctly.

end

