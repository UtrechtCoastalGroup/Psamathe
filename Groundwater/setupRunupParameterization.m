function [setup,runup,R2] = setupRunupParameterization(Hsig,Tp,beachSlope)

% SETUPRUNUPPARAMETERIZATION - Prediction of setup, runup and 2% exceedence
%                              value for natural beaches based on Stockdon
%                              et al. (2006, Coastal Engieering), page 587.
% Given time series of offshore significant wave height Hsig and wave
% period Tp, and a value for the beach slope, this function estimates time
% series of setup, runup, and the 2% waterlevel exceedence value based on
% the Stockdon et al. (2006) parameterizations, their Eqs.(9) to (12).
%
% INPUT
%   Hsig: offshore significant wave height, [Nt, 1]
%   Tp: offshore peak period, [Nt, 1]
%   beachSlope: the beach slope
% OUTPUT
%   setup, runup, R2: see above, [Nt, 1] each
%
% Gerben Ruessink, 14 May 2018, v1

% Some error checking
Hsig = Hsig(:);
Tp = Tp(:);
if length(Hsig) ~= length(Tp)
    error('Time series of wave height and period should be of equal length');
end

% offshore wave length and height
g = 9.81;                      % gravitational acceleration, m/s2
L0 = (g/(2*pi))*Tp.^2;         % deep-water wave length
H0 = Hsig;                     % deep-water wave height; we ignore the effect of de-shoaling / refraction

% apply parameterizations
setup = 0.35*beachSlope*sqrt(H0.*L0);                 % Eq. (10)
runup = sqrt(H0.*L0*(0.563*beachSlope^2 + 0.004))/2;  % Eq. (11) and (12)
R2 = 1.1*(setup + runup);                             % Eq. (9)

% check on dissipative conditions
IribarrenNumber = beachSlope ./ sqrt(H0./L0);
id = find(IribarrenNumber < 0.3);
R2(id) = 0.043*sqrt(H0(id).*L0(id));

% done!