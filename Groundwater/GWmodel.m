function [par,ts,sp,spt] = GWmodel(par,ts,sp)

% FUNCTION [par,ts,sp,spt] = GWmodel(par,ts,tp)
%
% The function GWmodel computes groundwater fluctuations in a sandy beach
% using the (non)linear Boussinesq equation - see Raubenheimer et al.
% (1999). More specifically, this function handles the time stepping and
% relies on BoussinesqGroundwater to perform the spatial central
% differences.

% INPUT
%    par, model parameters
%    ts, time series
%    sp, spatial (cross-shore) grids
% OUTPUT
%    as INPUT (with more fields)
%    spt, spatio-temporal information

% Create waterlevel vector on model time axis
if par.realtide == 1
    ts.tData = (ts.tAxis - ts.tAxis(1))*24*60*60;         % convert axis to seconds from first observation
    ts.tModel = 0:par.dt:ts.tData(end);                   % model time axis with time step (in s) of par.dt
    ts.tRealSU = interp1(ts.tData,ts.zShore,ts.tModel);   % shoreline elevation on the model time axis
    ts.tRealRU = interp1(ts.tData,ts.zRunup,ts.tModel);   % runup elevation on the model time axis
else
    error('Measured tide demanded. Abort.');
end

% Create bed levels on model grid
if par.realprofile == 1
    sp.x = 0:par.dx:sp.beach_location(end);                         % model domain 
    sp.profile = interp1(sp.beach_location,sp.beach_profile,sp.x);  % bed elevation
else
    error('Measured bed profile demanded. Abort.');    
end

% initialize several matrices
sp.zeta = zeros(size(sp.x));                                                   % water level
spt.zetat = zeros(ceil(ts.tModel(end)/par.outputtimes)+1,length(sp.zeta));     % output of water level
spt.tAxis = zeros(ceil(ts.tModel(end)/par.outputtimes)+1,1);                   % time axis
spt.shoreline = zeros(ceil(ts.tModel(end)/par.outputtimes)+1,1);               % time series of shoreline
spt.outcrop = NaN(ceil(ts.tModel(end)/par.outputtimes)+1,1);                   % time series of outcrop point
spt.x = sp.x;                                                                  % add x-axis so both time and spatial information is included in spt

% set water table for first instance
par.n = find(sp.profile > ts.tRealSU(1),1,'first');              % find shoreline location
sp.zeta(1:par.n-1) = ts.tRealSU(1);                              % up to shoreline, set to tidal level + setup  
sp.zeta(par.n:end) = ts.tRealSU(1):(par.onshorehead-ts.tRealSU(1))/(length(sp.zeta)-par.n):par.onshorehead; % linear increase from ts.tRealSU to par.onshorehead
sp.zeta(end) = sp.zeta(end-1);                                 % invoke d zeta/dx = zero at landward boundary 

q = waitbar(0,'Please wait');

% loop for all time steps
for i = 1:length(ts.tRealSU)

    waitbar(i/length(ts.tRealSU),q);
    
    par.zSL = ts.tRealSU(i);
    par.n = find(sp.profile > par.zSL,1,'first');                         % find shoreline location
    par.m = find(sp.profile < sp.zeta & sp.profile > par.zSL,1,'last');   % find outcrop point
    sp.tabledepth = sp.profile-sp.zeta;                                   % watertable depth  
    par.d = find(sp.tabledepth > par.minDepth,1,'first');                 % find first location exceeding minDepth (relevant to runup)

    % set boundary condition
    sp.zeta(1:par.n-1) = par.zSL;                            % up to shoreline, set to tidal level + setup
    sp.zeta(par.n:par.m-1) = sp.profile(par.n:par.m-1);      % set region between shoreline and outcrop point to bed height
    sp.zeta(end) = sp.zeta(end-1);                           % invoke d zeta/dx = zero at landward boundary 

    % use 4th order runge kutta method for time stepping
    z1 = BoussinesqGroundwater(sp.zeta,sp.profile,par);
    zeta1 = sp.zeta + 0.5.*z1.*par.dt;
    z2 = BoussinesqGroundwater(zeta1,sp.profile,par);
    zeta2 = sp.zeta + 0.5.*z2.*par.dt;
    z3 = BoussinesqGroundwater(zeta2,sp.profile,par);
    zeta3 = sp.zeta + z3.*par.dt;
    z4 = BoussinesqGroundwater(zeta3,sp.profile,par);
    sp.zeta = sp.zeta + 1/6*(z1+2*z2+2*z3+z4).*par.dt;

    % estimate run-up infiltration
    if par.runup
       par.zRU = ts.tRealRU(i);
       par.r = find(sp.profile > par.zRU,1,'first');                                     % find runup location
       sp.infiltrationFunction = (sp.x - sp.x(par.d))/(sp.x(par.r) - sp.x(par.d));       % dimensionless function
       sp.infiltrationFunction(1:par.d) = 0;                                             % which should be zero up to location of minDepth
       sp.infiltrationFunction(par.r:end) = 0;                                           % and above the runup location  
       sp.Ul = par.Cl*(par.K/par.ne)*sp.infiltrationFunction;                            % infiltration velocity averaged over numerous uprush/backwash cycles (= infiltration flow rate per unit area)  
       sp.zeta = sp.zeta + sp.Ul.*par.dt;                                                % add it to sp.zeta  
    end;

    % export data
    spt.zetat(ceil(ts.tModel(i)/par.outputtimes)+1,:) = sp.zeta;                            % water table
    spt.tAxis(ceil(ts.tModel(i)/par.outputtimes)+1) = ts.tModel(i)/24/60/60+ts.tAxis(1);    % in Matlab time
    spt.shoreline(ceil(ts.tModel(i)/par.outputtimes)+1) = sp.x(par.n);                      % time series of shoreline position
    if ~isempty(par.m)
       spt.outcrop(ceil(ts.tModel(i)/par.outputtimes)+1) = sp.x(par.m);                     % time series of outcrop point
    end;

end

close(q);

% done!