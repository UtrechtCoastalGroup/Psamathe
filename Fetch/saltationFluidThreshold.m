function uStar_it = saltationFluidThreshold(par)
%SALTATIONFLUIDTHRESHOLD Computation of saltation fluid treshold according
%to Shao and Lu (2000)
%
% This function file computes the saltation fluid threshold as proposed by
% Shao and Lu in their 2000 JGR paper. It is a simplification to the
% equations proposed by Greeley and Iverson (1985), which is an extension
% to the traditional Bagnold threshold equation. For sand this extension is
% not that relevant, by the way. The code assumes dry sand.
%
% INPUT
%   par, a structure with parameter values. Of relevance here are the
%        fields AN, rhoS, rhoA, g, D50, gamma, beachSlope and angleOfRepose.
%        Typical values are AN = 0.111, rhoS = 2650 kg/m^3, rhoA = 1.25 kg/m^3,
%        g = 9.81 m/s^2,gamma = 2.9 × 10^-4 Nm^-1 and angleOfRepose = 33 deg.
%        This gamma is based on De Kok et al. (2012). D50 is the median grainsize (m).
% OUTPUT
%   uStar_it, the saltationFluidTreshold
%
% Example: With the above values and D50 = 250e-6 m, uStar_it = 0.278 m/s on a horizontal bed.
%
% REFERENCES
% Shao, Y. and Lu, H., 2000. A simple expression for wind erosion threshold
%   friction velocity. Journal of Geophysical Research, 105, 22,437-22,443.
% De Kok, J.F., E.J.R. Partelli, T.I. Michaels and D.B. Karam, 2012. The
%   physics of wind-blown sand and dust. Rep. Prog. Phys, 75, 72 pp.
%
% v1.0, Gerben Ruessink, April 1, 2019
% v1.1, Gerben Ruessink, October 7, 2019. Added bed slope effect according
%                        to Iversen and Rasmussen (1994)

% error checking
fieldPresence = isfield(par,{'AN','rhoS','rhoA','g','D50','gamma','beachSlope','angleOfRepose'});
if sum(fieldPresence) ~= 8
    error('One or more parameters missing in saltationFluidTreshold.m. Required: AN, rhoS, rhoA, g, D50, gamma, beachSlope,angleOfRepose.');
end

% bedslope effect according to Iversen and Rasmussen (1994), see also Van
% Dijk et al. (1999), Eq. (10)
bedslopeEffect = sqrt(cos(par.beachSlope) + (sin(par.beachSlope)/tand(par.angleOfRepose)));

% compute saltation fluid threshold according to Shao and Lu (2000), see
% also Eq. (2.8) in De Kok et al. (2012)
uStar_it = bedslopeEffect * par.AN * sqrt(((par.rhoS - par.rhoA)/par.rhoA)*par.g*par.D50 + ...
        par.gamma/(par.rhoA*par.D50));

end

