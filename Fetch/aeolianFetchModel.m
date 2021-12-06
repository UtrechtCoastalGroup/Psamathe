function [ts,spt] = aeolianFetchModel(ts,spt,sp,par)

%AEOLIANFETCHMODEL AEOLIAN TRANSPORT ACROSS A BEACH (version 1.6)
%   In 2002, Bauer and Davidson-Arnott published their now classic, conceptual
%   work on (supply-limited) aeolian transport as a function of wind angle, beach
%   geometry and fetch effects. Later on, Delgado-Fernandez tried to
%   implement the conceptual model in an actual set of equations. Here, we
%   build on her work and extend it as follows:
%   (1) Surface moisture is not a priori measured, but predicted from a
%   ground water - retention curve model.
%   (2) Surface moisture can vary spatially.
%   The latter point implies that the cross-shore evolution of the aeolian
%   transport rate is no longer a (simple) function, but rewritten in
%   advection - pickup form.
%
%   We assume here that:
%   (1) no aeolian transport for alongshore or offshore directions.
%   (2) no aeolian transport if surface moisture content exceeds a certain
%   threshold, typically 10% (= input parameter)
%
%   If the wind speed is large enough and the direction is (obliquely)
%   onshore, then
%   q(x) = 0, if surfMoist > par.moistMax
%   q(x) = min(qPotential,q(x-1) + qPotential*sin((pi/2)*(F/Fc(x))).
%
%   -- qPotential is the potential transport, computed here according to Hsu
%   (1971), Kok et al. (2012) or Lettau and Lettau (1978).
%   -- F is the downwind distance ('fetch') over which critical fetch is
%   constant. So every the critical fetch changes, the 'fetch' is reset.
%   -- Fc is the critical fetch, which depends on wind speed and moisture
%   content. See comments with the internal function.
%
%   As input the code requires:
%   -- spt, spatial-temporal matrices. The relevant one here is
%   spt.surfMoist, that is, cross-shore surface moisture profiles for a
%   given set of time steps.
%   -- ts, temporal series. The relevant ones here are wSpeed, wDirBeach
%   and wDirDunefoot. The directions are the wind direction on the beach
%   and at the dunefoot. We assume that the wind, when it picks up the
%   sand, has the same direction as the 'regional' wind, but that it has
%   steered more alongshore ('wDirDunefoot') at the dunefoot. The latter is
%   relevant in the so-called cosine effect. We assume (not checked!) that
%   the time series of wind charcteristics are at the same time axis as
%   spt.surfMoist. This is thus a preprocessing step! (i.e., not done here
%   in the code). Note that wDirBeach and wDirDunefoot should be in degrees
%   relative to the shore normal, so angles that allow for aeolian
%   transport are between (but not including) -90 and +90 degrees.
%   -- sp, spatial matrices. Relevant here are x (cross-shore grid) and
%   profile (= bed profile). Used to get the cross-shore grid size and to
%   determine for which part of the profile aeolian transport needs to be
%   computed at all (see par.zUp below).
%   -- par, parameters. Relevant here are aeolianModel ('Hsu', 'Kok' or 
%   'Lettau'), moistMax, rainIntensityMax and zUp. A considerable number of
%   other parameters are used in functions called from this function. zUp
%   is the elevation above which we assume that sedimentation is of aeolian origin.
%
%   The output provided is:
%   -- within ts, the matrices qPotential (in kg/m/s),
%   qPotentialCosine being equal to qPotential but the part that may go
%   from the beach on to the dune (cosine effect, using wDirDunefoot), and qActual
%   and qActualCosine being the predicted actual transport on to the dune
%   (the latter including the cosine effect with wDirDunefoot).
%   -- within spt, the matrices F and Fc (in m) being the local and critical fetch,
%   respectively, and qCum (in km/m/s) being the cross-shore evolution of the aeolian transport.
%   The final value of qCum times cos(wDirDunefoot) equals qActualCosine in ts.
%
%   Some more remarks:
%   -- in contrast to earlier versions, the present version does not
%   contain a minimum wind speed. In the Kok et al. (2012) and Lattau
%   and Lattau (1978) models an initiation of motion is included. In Hsu
%   (1971) this is not the case, but can be 'activated' with
%   par.thresholdWind = 1.
%
% v1.0, Gerben Ruessink, February 2018
% v1.1, Gerben Ruessink, August 2018 - added qActual (without cosine)
% v1.2, Gerben Ruessink, September 2018 - corrected bug when all moisture
%       contents are below par.moistMax.
% v1.3, Gerben Ruessink, December 2018 - modified Fc slightly. In
%       hindsight, the coding in the fetch part was incorrect because of the
%       Fc modifications. This version should not be used.
% v1.4, Gerben Ruessink, April 2019 - extended with multiple models for
%       potential transport and removed the wMinSpeed criterium. Also set
%       Fc equal to Delgado-Fernandez formulation, and hence modified
%       the fetch part. The rounding of moisture to multiples of 0.5 has
%       disappeared. Included error checking on par. Recoded v1.2
%       implementation for efficiency. Fetch is now an output in spt too.
% v1.5, Gerben Ruessink, May 2019 - extended to include rain in two ways.
%       First, if it rains, then the critical uStar is raised by 35%. This
%       is not done in this code, but in the code that computes the potential
%       transport rates. Second, if the intensity exceeds a given maximum,
%       qActual is set to 0.
% v1.6, Gerben Ruessink, December 2021 - corrected various typos.

% error checking
fieldPresence = isfield(par,{'aeolianModel','moistMax','zUp','rainIntensityMax'});
if sum(fieldPresence) ~= 4
    error('One or more parameters missing in aeolianFetchModel.m. Required: aeolianModel, moistMax, rainIntensityMax and zUp.');
end

% Initialize output matrices, to be stored in ts
qPotential = zeros(size(ts.wSpeed));
qPotentialCosine = zeros(size(ts.wSpeed));
qActual = zeros(size(ts.wSpeed));
qActualCosine = zeros(size(ts.wSpeed));

% Initialize output matrices, to be stored in spt
Fc = NaN(size(spt.surfMoist));
fetch = Fc;
qCum = zeros(size(spt.surfMoist));

% The code assumes that x increases in the shoreward direction. Let's find
% the most landward point up to which computations need to be performed
% based on par.zUp. Typically, par.zUp will be between 2 or 3 m +MSL. We
% also assume here that the profile crosses par.zUp only once.
id = find(sp.profile < par.zUp);
if length(id) == length(sp.profile) || isempty(id) || length(id) < 2
    error('aeaolianFetchModel:zUp','Error: par.zUp is not found on the beach profile.');
end
nUp = id(end);

% cross-shore grid size
dx = abs(sp.x(2)-sp.x(1));

% number of time steps
Nt = length(ts.wSpeed);

% work through all time steps
for i = 1:Nt
   
    % check on wind direction; should have onshore component
    if abs(ts.wDirBeach(i)) >= 90
       continue;
    end
    
    % is there any rain?
    if ts.rainIntensity(i) > 0
        par.rainPotential = 1;     % yes, include effect of rain on qPotential.
    else
        par.rainPotential = 0;     % no, it is dry.
    end
    
    % potential transport and its cross-shore component
    switch par.aeolianModel
        case 'Hsu'
            qPotential(i) = aeolianTransportRateHsu(ts.wSpeed(i),par);
        case 'Kok'
            qPotential(i) = aeolianTransportRateKok(ts.wSpeed(i),par);
        case 'Lettau'
            qPotential(i) = aeolianTransportRateLettau(ts.wSpeed(i),par);
        otherwise
            error('Unknown aeolian transport model');
    end
    qPotentialCosine(i) = qPotential(i).*cosd(ts.wDirForedune(i));

    % No need to compute anything if wind speed was too low
    if qPotential(i) == 0
        continue;
    end
    
    % Now we can proceed with actual transport computations. Well, we need
    % to check first whether there is sufficient dry sand. If this fails, it
    % represents the extreme case when qPotential is non-zero, but qActual
    % is zero.
    surfMoist = spt.surfMoist(i,1:nUp);
    if all(surfMoist >= par.moistMax) || sum(isnan(surfMoist)) == nUp
        continue;
    end
    
    % And we need to check for rain
    if ts.rainIntensity(i) >= par.rainIntensityMax
        continue;
    end
        
    % We have reached the stage where we can compute actual transport
    % rates. Let's first set all points seaward of most landward
    % par.moistMax to NaN.
    idTooMoist = find(isnan(surfMoist) | surfMoist > par.moistMax);
    % In rare situations it may happen that all surfMoist are less than
    % par.moistMax. Then, without the isnan(surfMoist) in the previous line
    % idTooMoist would be an empty matrix and the code below would crash.
    surfMoist(idTooMoist) = NaN;
    
    % Compute critical fetch for this situation
    Fc(i,1:nUp) = criticalFetch(surfMoist,ts.wSpeed(i));
    
    % downwind grid size
    dxWind = dx / cosd(ts.wDirBeach(i));
    
    % We will use the unique function to get the parts with constant
    % critical fetch. Note that each NaN is considered unique.
    % Hence we start the loop s at idTooMoist(end)+1, as this is the first
    % unique non-NaN value.
    [uniqueValues,~,IC] = unique(Fc(i,1:nUp),'stable');
    for s = idTooMoist(end)+1:length(uniqueValues)
       
        % These points have equal surface moisture
        id = find(IC==s);
        
        % Then this is the fetch F. It is 1 entry longer than id. We assume
        % now that 0 is at the previous (upwind) grid point. Hence, in the
        % computation of qCum below we use F(2:end).
        F = 0:dxWind:length(id)*dxWind;
        
        % Limit F to Fc; otherwise the sin-term below gets below 1.
        F(F > Fc(i,id(1))) = Fc(i,id(1));
        
        % Compute qCum for this stretch. It is q(x-1) +
        % qPot*sin((pi/2)*(F/Fc)), up to a maximum of qPot
        qCum(i,id) = min(qPotential(i),qCum(i,id(1)-1) + qPotential(i)*sin((pi/2)*(F(2:end)./Fc(i,id))));
        
        % Store fetch too (in hindsight)
        fetch(i,id) = F(2:end);
        
    end
    qCum(i,nUp+1:end) = NaN; % Set everything above nUp to NaN, is better than 0.

    % Store actual amount and then apply cosine effect for z = zUp.
    qActual(i) = qCum(i,nUp);
    qActualCosine(i) = qCum(i,nUp).*cosd(ts.wDirForedune(i));
        
end

% store all output in either ts or spt
ts.qPotential = qPotential;
ts.qPotentialCosine = qPotentialCosine;
ts.qActual = qActual;
ts.qActualCosine = qActualCosine;
spt.Fc = Fc;
spt.qCum = qCum;
spt.F = fetch;

end

function Fc = criticalFetch(surfMoist,wSpeed)
%CRITICALFETCH Compute critical fetch Fc
%   We compute the critical fetch Fc as:
%   Fc = alpha * [4.38*U - 8.23], where U is the wind speed wSpeed and
%   alpha is a correction factor that depends on the surface moisture
%   content surfMoist as specified in Delgado-Fernandez (2011), with a
%   modification for moisture contents in excess of 10%.

% alpha
alpha = ones(size(surfMoist));
alpha(surfMoist >= 4 & surfMoist < 6) = 1.25;
alpha(surfMoist >= 6 & surfMoist <= 10) = 1.75;
alpha(surfMoist > 10) = 2.50;
alpha(isnan(surfMoist)) = NaN;

% critical fetch
Fc = alpha.*(4.38*wSpeed - 8.23);

end

