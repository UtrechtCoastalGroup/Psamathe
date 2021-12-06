function qLettau = aeolianTransportRateLettau(U,par)
%AEOLIANTRANSPORTRATELettau Lettau-Lettau equation for aeolian transport rate
%   This function computes the (dry-sand) potential aeolian transport rate
%   for a given velocity and a parameter that links the wind velocity to a
%   shear velocity. The computation follows Lettau and Lettau (1978) as 
%   described in Sherman et al. (2013), and uses the Shao and Lu (2000) 
%   approach for the saltation fluid threshold.
%
% INPUT
%   U, wind velocity, N values (row or column) of wind velocity (m/s)
%   par, a structure with parameter values. Of relevance here are the
%      fields CL, rhoA, g, D50 and a. Typical values are CL = 6.7, rhoA =
%      1.25 kg/m^3, g = 9.81 m/s^2 and a = 0.04. The 'a' relates the wind
%      velocity at a specific height to the shear velocity assuming a
%      logarhithmic velocity profile. See papers by Hsu on this topic.
%      Consider Sherman et al. (2013) for a 'calibrated' CL.
%
% OUTPUT
%   qLettau, steady-state (potential) aeolian transport rate, kg/m/s for
%      each U.
%
% REFERENCES
% Lettau and Lettau (1978) in Sherman et al. (2013). I don't have the
% original Lettau and Lettau book chapter.
% Sherman, D.J. et al., 2013. Recalibrating aeolian sand transport models.
%    Earth Surface Processes and Landforms, 38, 169-178.
%
% v1.0, Gerben Ruessink, April 2, 2019
% v1.1, Gerben Ruessink, May 14, 2019. Included rain effect (in uStar_it).

% error checking
fieldPresence = isfield(par,{'CL','rhoA','g','a','D50','rainPotential','CRain'});
if sum(fieldPresence) ~= 7
    error('One or more parameters missing in saltationFluidTreshold.m. Required: CL, rhoA, g, a, D50, rainPotential and CRain.');
end

% saltation fluid trenshold based on the method of Shao and Lu (2000)
uStar_it = saltationFluidThreshold(par);

% include effect of rain if so desired based on Arens (1996)
if par.rainPotential
   uStar_it = uStar_it*(1+par.CRain);    
end

% shear velocity
uStar = par.a*U;

% steady-state transport rate, Eq. (2.34) in De Kok et al. (2012). If
% uStar_it > uStar, qLettau will become less than 0. This is, of course,
% prevented.
qLettau = max(0,par.CL*sqrt((par.D50/250e-6))*(par.rhoA/par.g)*uStar.^2.*(uStar-uStar_it));

end

