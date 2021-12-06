function qKok = aeolianTransportRateKok(U,par)
%AEOLIANTRANSPORTRATEKOK Kok equation for aeolian transport rate
%   This function computes the (dry-sand) potential aeolian transport rate
%   for a given velocity and a parameter that links the wind velocity to a
%   shear velocity. The computation follows Kok et al. (2012), their Eq.
%   (2.34), and uses the Shao and Lu (2000) approach for the saltation
%   fluid threshold.
%
% INPUT
%   U, wind velocity, N values (row or column) of wind velocity (m/s)
%   par, a structure with parameter values. Of relevance here are the
%      fields CDK, rhoA, g and a. Typical values are CDK = 5,  rhoA =
%      1.25 kg/m^3, g = 9.81 m/s^2 and a = 0.04. The 'a' relates the wind
%      velocity at a specific height to the shear velocity assuming a
%      logarhithmic velocity profile. See papers by Hsu on this topic.
%
% OUTPUT
%   qKok, steady-state (potential) aeolian transport rate, kg/m/s for
%      each U.
%
% REFERENCES
% Kok, J.F., E.J.R. Partelli, T.I. Michaels and D.B. Karam, 2012. The
%   physics of wind-blown sand and dust. Rep. Prog. Phys, 75, 72 pp.
%
% v1.0, Gerben Ruessink, April 1, 2019
% v1.1, Gerben Ruessink, May 13, 2019. Included rain effect (in uStar_it).
% v1.2, Gerben Ruessink, December 6, 2021. Corrected several typos.

% error checking
fieldPresence = isfield(par,{'CDK','rhoA','g','a','rainPotential','CRain'});
if sum(fieldPresence) ~= 6
    error('One or more parameters missing in aeolianTransportRateKok.m. Required: CDK, rhoA, g, a, rainPotential and CRain.');
end

% saltation fluid trenshold based on the method of Shao and Lu (2000)
uStar_it = saltationFluidThreshold(par);

% include effect of rain if so desired based on Arens (1996)
if par.rainPotential
   uStar_it = uStar_it*(1+par.CRain);    
end

% shear velocity
uStar = par.a*U;

% steady-state transport rate, Eq. (2.34) in Kok et al. (2012). If
% uStar_it > uStar, qKok will become less than 0. This is, of course,
% prevented.
qKok = max(0,(par.CDK*par.rhoA/par.g)*uStar_it*(uStar.^2 - uStar_it^2));

end
