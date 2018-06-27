% This script file provides an example how to run the groundwater - surface
% moisture model given a measured bed profile and hydrodynamical forcing
% (tides, waves) at the seaward boundary. The example is taken from the
% MegaPex experiment on the Sand Engine in the Netherlands, and uses
% optimum parameter values given (previous) data-model optimization. The
% example is presented and analysed in:
%
% Brakenhoff, L.B., Y. Smit, J.J.A. Donker and G. Ruessink, under review.
% Tide-induced variability in beach surface moisture: observations and
% modelling. Earth Surface Processes and Landforms.
%
% Hopefully there are sufficient comments to guide you through the code. As
% a start, all input-output is stored into four major structures:
% - sp, for spatial information (independent of time)
% - ts, for time series at a given location
% - spt, for spatio-temporal matrices
% - par, for all parameter values and model choices
%
% Contact person: Gerben Ruessink, b.g.ruessink@uu.nl
% Feel free to use in your own project. Please refer to the paper above in
% any publication that may arise from your work.
%
% v1, Gerben Ruessink, 27-06-2018; works on MATLAB Version: 9.3.0.713579 (R2017b)

clear
close all
clc

% Add path where the code is. Here I assume that the example is run from an
% example subdirectory within the path with the code (as on GitHub).
examplePath = strrep(which('runGWLMoistureModel'),'runGWLMoistureModel.m', '');
codePath = strrep(examplePath,['example' filesep],'');
addpath(codePath);

%% (1) input for beach profile and waterlevels (tide, setup, R2)

% Let's first make the time-axis for the model input - this is later on
% interpolated on the actual computation axis
startTime = datenum([2014 9 16 0 0 0]);                   % This is the approximate start of the MegaPex campaign. 
endTime = datenum([2014 10 21 12 0 0]);                   % This approximately corresponds to the last low-tide with good data.
tStep = 10;                                               % Time step (in minutes)
tAxis = (startTime:datenum([0 0 0 0 tStep 0]):endTime)';  % Time vector for waterlevels
ts.tAxis = tAxis;                                         % We store all time series in structure ts

% The beach profile: x and z, with x = 0 at the offshore boundary and
% positive onshore. See header of the file that is loaded in the next line.
xz = load([examplePath,'bedProfiles.txt']);
par.realprofile = 1;
sp.beach_location = xz(:,1);
sp.beach_profile = xz(:,end);

% The shoreline elevation will be based on zTide and setup. We use the 'offshore'
% tidal values and the estimated setup. The runup elevation equals R2. We
% use the parameterizations of Stockdon et al. (2006). 'offshore' implies
% the data at the PT.
par.realtide = 1;
zetaOffshore = load([examplePath,'zetaOffshore.txt']);
tZetaOffshore = datenum(zetaOffshore(:,1:6));
zetaOffshore = interp1(tZetaOffshore,zetaOffshore(:,end),ts.tAxis);

wavesOffshore = load([examplePath,'wavesOffshore.txt']);
tWavesOffshore = datenum(wavesOffshore(:,1:6));
H0 = interp1(tWavesOffshore,wavesOffshore(:,7),ts.tAxis);
T0 = interp1(tWavesOffshore,wavesOffshore(:,8),ts.tAxis);

p = polyfit(sp.beach_location,sp.beach_profile,1);
beachSlope = p(1);
[setup,runup,R2] = setupRunupParameterization(H0,T0,beachSlope);

% zShore = offshore waterlevel + setup
ts.zShore = zetaOffshore + setup;

% zRunup = offshore waterlevel + setup + runup. For (setup + runup) we use R2.
ts.zRunup = zetaOffshore + R2;

%% (2) Groundwater model settings
par.dt = 2;                    % Time step in seconds
par.dx = 0.5;                  % Grid size in meters
par.nl = 1;                    % Include non-linear parts in groundwater model (1 = Yes; 0 = No)
par.runup = 1;                 % Include runup infiltration flux (1 = Yes; 0 = No)
par.onshorehead = 1;           % Initial overheight onshore
par.outputtimes = tStep*60;    % Generate output every # seconds
par.D = 7;                     % Aquifer depth
par.K = 7.7868e-04;            % Hydraulic conductivity [optimum value based on data-model optimization]
par.ne = 0.3;                  % Effective prosity
par.Cl = 0.2164;               % Infiltration coefficient [optimum value based on data-model optimization]
par.minDepth = 0.2;            % Minimum water table depth in runup-infiltration

%% (3) Run the groundwater model
[par,ts,sp,spt] = GWmodel(par,ts,sp);

% The relevant output is in the matrix spt:
% .zetat: groundwater levels. Note that if a location was submerged, the
%         value equals zShore.
% .tAxis: time axis
% .x: cross-shore axis

%% (4) Settings for soil-water-retention curve

% Parameters in Van Genuchten soil-water-retention curve [values based on
% Egmond data]
par.thetaSat = 25.1951;
par.thetaRes = 4.2215;
par.alpha= 5.3072;
par.n = 3.1756;
par.m = 1 - 1/par.n;

%% (5) Run the surface moisture model
spt = surfaceMoisture(par,sp,spt);

% The relevant output is in the matrix spt (obviously):
% .surfMoist: surface moisture content

%% (6) Some post-processing, save output and a make/save a figure

% Set all moisture values below xShorelineSetup to NaN
for i = 1:length(spt.shoreline)
    id = findnearest(sp.x,spt.shoreline(i));
    spt.surfMoist(i,1:id) = NaN;
end

save([examplePath,'exampleOutput.mat'],'par','ts','sp','spt');

% Plot spatio-temporal moisture content for October 11-20, 2014 (as in the Brakenhoff et al. paper). 

figure(1);  
subplot 211 
startDate = datenum([2014 10 11 0 0 0]);
endDate = datenum([2014 10 21 0 0 0]);
pcolor(spt.tAxis,sp.x,spt.surfMoist'); shading('flat');
set(gca,'xlim',[startDate endDate],'xtick',startDate:endDate,'xticklabel',11:21,...
    'ylim',[0 100],'ytick',0:20:100,'fontname','times','fontsize',11,'clim',[par.thetaRes par.thetaSat],...
    'tickdir','out','box','on');
colormap(flipud(colormap))
H = colorbar;
ylabel(H,'Moisture content (%)','fontsize',11,'fontname','times');
ylabel(gca,'Cross-shore distance (m)','fontsize', 11,'fontname','times');
xlabel(gca,'Time (day in October 2014)','fontsize', 11,'fontname','times');

print([examplePath,'exampleSurfaceMoistureOutput.png'],'-dpng','-r600');
