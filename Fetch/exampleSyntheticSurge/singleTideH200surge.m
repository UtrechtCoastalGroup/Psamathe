% Example script file to apply the aeolian fetch model to a synthetic case
% of a 1:30 beach and a 2-m 12-hour tide with a 2-m surge for a
% shore-normal wind with a speed of 17.5 m/s.
%
% Gerben Ruessink, December 6, 2021

clear
close all
clc

% add paths - adjust!
addpath(genpath('..\Psamathe-master\'));

%% Time, beach profile and water levels

% Implement a 1:30 slope, starting at x = 0 and z = -2 m, with all values
% above 5 m set to 5 m. The profile is sufficiently long at the landward
% side.
par.realprofile = 1;                         
par.beachSlope = 1/30;
x = 0:400;
z = min(x*par.beachSlope-2,5);
sp.beach_location = x;
sp.beach_profile = z;

% Make a 12-hour tide with a range of 2 m. This tide is repeated
% 50 times (= 25 days), with a 10-min resolution. Then, we make a 'fake'
% Matlab time axis and add the surge.
par.realtide = 1;
T = 12/24; 
H = 2;
dt = (10/60)/24;
Ndays = 25;
t = 0:dt:Ndays-dt;
eta = ((H/2)*cos(2*pi*t/T))';
startTime = datenum([1 1 1 0 0 0]);
endTime = datenum([1 1 Ndays 23 50 0]);
tAxis = (startTime:dt:endTime)';   
ts.tAxis = tAxis;
ts.rainIntensity = zeros(length(tAxis),1);  % no rain
sigma = 0.5;
mu = 19.5;
amplitude = 2;
surgelevel = amplitude*exp(-0.5*(((ts.tAxis-ts.tAxis(1))-mu)/sigma).^2);
ts.zShore = eta + surgelevel;
ts.zRunup = ts.zShore;       % runup will be switched off in model settings

%% Model settings

% No wave runup, values loosely based on Egmond case study
par.dt = 5;
par.dx = 0.5;
par.nl = 1;
par.runup = 0;
par.onshorehead = 0.5;
par.outputtimes = 600;
par.K = 40 / (3600*24);
par.ne = 0.3;
par.D = 15;
par.Cl = 0.5;                  
par.minDepth = 0.2;            

%%  Run ground water module
[par,ts,sp,spt] = GWmodel(par,ts,sp);

% The relevant output is in the matrix spt:
% .zetat: groundwater levels. Note that if a location was submerged, the
%         value equals zShore.
% .tAxis: time axis
% .x: cross-shore axis

%% Run surface moisture module. Base alpha and n based on Tuller and Or.
par.thetaSat = 20;
par.thetaRes = 2;
par.alpha = 3.5;
par.n = 3.2;
par.m = 1 - 1/par.n;
spt = surfaceMoisture(par,sp,spt);
for i = 1:length(spt.shoreline) % Set all moisture values below xShorelineSetup to NaN
    id = findnearest(sp.x,spt.shoreline(i));
    spt.surfMoist(i,1:id) = NaN;
end

% The relevant output is in the matrix spt (obviously):
% .surfMoist: surface moisture content

%% Run fetch model with Kok et al. (2012)

% set all relevant parameters
par.moistMax = 10;
par.zUp = 2.5;
par.D50 = 250e-6;
par.a = 0.04;
par.g = 9.81;
par.aeolianModel = 'Kok';
par.AN = 0.1109;
par.CDK = 5;
par.rhoA = 1.25;
par.rhoS = 2650;
par.gamma = 2.9e-4;
par.angleOfRepose = 33;
par.rainIntensityMax = 1000;  % not relevant
par.CRain = 0;                % not relevant

% make wind axes
U = 17.5;
DirSN = 0;
ts.wSpeed = U*ones(length(ts.tAxis),1);
ts.wDirBeach = DirSN*ones(length(ts.tAxis),1);
ts.wDirForedune = ts.wDirBeach;

% and run the model
[ts,spt] = aeolianFetchModel(ts,spt,sp,par);
        
% The relevant output is in the matrix ts:
% .qPotential: potential transport rate (in kg/m/s)
% .qActual: actual transport rate at the beach-dune intersection (in
% kg/m/s)
% NOTE: wind direction was 0 here (shore-normal winds); so qPotentialCosine
% and qPotentialActual are the same as qPotential and qActual.
