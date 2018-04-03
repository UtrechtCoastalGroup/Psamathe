function ddt = BoussinesqGroundwater(eta,zBed,par)

% FUNCTION DDT = BOUSSINESQGROUNDWATER(ETA,ZBED,PAR)
% 
% BoussinesqGroundwater contains the Boussinesq equation for
% finite-amplitude (nonlinear) water table fluctuations for 
% constant aquifer thickness, Eq. (1) from Raubenheimer et al. (1999)
%
% d(eta)/dt = (K/ne) d{[D+eta]d(eta)/dx)}dx
%
% This can be rewritten as:
%
% d(eta)/dt = (KD/ne)d2eta/dx2 + (K/ne) d{eta d(eta)/dx} dx
%
% For small-amplitude (linear) water table fluctuations, the second term on
% the right-hand side can be ignored.
%
% INPUT
%   eta, water table height (x) [m]
%   zBed, bed elevation (x) [m]
%   par contains the following fields:
%      .zSL, shoreline water level [m]
%      .K, hydraulic conductivity in the saturated beach
%      .D, aquifer thickness [m]
%      .ne, effective porosity
%      .dx, (constant) grid size [m]
%      .nl, boolean to include non-linear term (1) or not (0)
% OUTPUT
%   ddt, temporal gradient in eta at all x

n = find(zBed > par.zSL,1,'first');                  % find shoreline
m = find(zBed < eta & zBed > par.zSL,1,'last');      % outcrop point
eta(1:n-1) = par.zSL;                                % set boundary at shoreline
eta(n:m-1) = zBed(n:m-1);                            % set groundwater level to fully wet between shoreline and outcrop point
eta(end) = eta(end-1);                               % invoke d eta/dx=zero at the edge

% initialize various terms
ddx = zeros(size(zBed));
ddx2 = zeros(size(zBed));
ddt = zeros(size(zBed));

% second difference operator in space (= d2eta/dx2 in first term rhs)
ddx(2:end-1) = ((eta(3:end)-2*eta(2:end-1)+eta(1:end-2))) / par.dx^2;

% eta d(eta)/dx (central difference)
ddx2(2:end-1) = eta(2:end-1).*((eta(3:end)-eta(1:end-2)) / (par.dx));
% d{}/dx of previous line (= d{eta d(eta)/dx} dx in second term rhs)
ddx2(2:end-1) = ((ddx2(3:end)-ddx2(1:end-2)) / (par.dx));

% combine the two, multiply with nonlin (if 0, then linear)
ddt(2:end-1) = par.K.*par.D/par.ne*ddx(2:end-1) + par.K./par.ne*ddx2(2:end-1).*par.nl;

end