function spt = surfaceMoisture(par,sp,spt)

% FUNCTION SPT = SURFACEMOISTURE(PAR,SP,SPT)
%
% The function surfaceMoisture estimates the surface moisture content (by mass)
% using the Van Genuchten soil water retention curve, that is, Equations (21) and
% (22) in Van Genuchten (1980), Soil Sci. Soc. Am. J., 44, 892-898.
%
% INPUT
%   par, structure with model parameters. Relevant here are:
%        thetaRes = residual moisture content (in %)
%        thetaSat = saturated moisture content (in %)
%        alpha = fit parameter
%        n = fit parameter
%        m = 1 - 1/n
%  sp, structure with spatial information. Of relevance here is the
%        cross-shore profile, sp.profile
%  spt, structure with spatial-temporal information. Of relevance here are
%        the predicted groundwater levels, spt.zetat. This is output from
%        GWmodel.m
%
% OUTPUT
%  spt, as input, but with surfMoist added.
%
% Gerben Ruessink, v1, March 22, 2017. Based on work by Jasper Donker.

% Compute depth beneath bed profile for all time steps simulteneously
Nt = size(spt.zetat,1);
hBelowBed = ones(Nt,1)*sp.profile - spt.zetat;

% Values below 0 implies inundated (by the tide): set these to 0.
hBelowBed(hBelowBed < 0) = 0;

% Apply Van Genuchten (1980) curve given parameters in par structure. Note
% that hBelowBed == 0 results in surfMoist = par.thetaSat;
spt.surfMoist = par.thetaRes + (par.thetaSat - par.thetaRes) ./ ...
    (1 + (par.alpha*hBelowBed).^par.n).^par.m;

% done!

