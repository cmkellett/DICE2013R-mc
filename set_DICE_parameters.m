function Params = set_DICE_parameters(N)
%
% Params = set_DICE_parameters(N)
%
% N = horizon length of problem
%
% Function to set parameters and exogenous signals
%
% -------------------------------------------------------------------------
%  This file is part of DICE2013R-mc.
% 
% DICE2013R-mc -- Dynamic Integrated model of Climate and Economy 
%                   2013R - Matlab and Casadi.
% Copyright (C) 2016 Timm Faulwasser, Christopher M. Kellett, and 
%       Steven R. Weller. 
% 
% DICE2013R-mc can be downloaded from https://github.com/cmkellett/DICE2013R-mc.
% 
% DICE2013R-mc is free software; you can redistribute it and/or modify it 
% under the terms of the GNU Lesser General Public License (LGPL) as 
% published by the Free Software Foundation; either version 3 of the 
% License, or (at your option) any later version.
% 
% DICE2013R-mc is distributed in the hope that it will be useful, but 
% without any warranty; without even the implied warranty of merchantability 
% or fitness for a particular purpose.  See the GNU Lesser General Public 
% License for more details.                        
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Horizon length
Params.N = N;
% -------------------------------------------------------------------------
% Exogenous Signals

Params.sigma = zeros(1,N+1);
Params.L = zeros(1,N+1);
Params.A_TFP = zeros(1,N+1);
Params.E_Land = zeros(1,N+1);
Params.F_EX = zeros(1,N+1);
Params.theta1 = zeros(1,N+1);

sigma(1)    = 0.549128362880230; % Initial sigma (GAMS sig0)
L(1)        = 6838;           % Initial population in millions (GAMS pop0)
A_TFP(1)    = 3.8;            % Initial TFP (GAMS a0)

for i=1:N+1
    sigma(i+1) = sigma(i) * exp(-0.01 * ((0.999^5)^i) * 5);
    L(i+1) = L(i) * (10500/L(i))^0.134;
    A_TFP(i+1) = A_TFP(i) / (1 - 0.079 * exp(-0.006 * 5 * (i-1)));
    E_Land(i) = 3.3*(0.8^(i-1));
    if (i < 19)
        F_EX(i) = 0.25 + 0.025*(i-1);
    else
        F_EX(i) = 0.25+0.45;
    end
    theta1(i) = sigma(i) * (344/2800) * 0.975^(i-1);
end

Params.sigma = sigma;
Params.L = L;
Params.A_TFP = A_TFP;
Params.E_Land =E_Land;
Params.F_EX = F_EX;
Params.theta1 = theta1;

% -------------------------------------------------------------------------
% Constants
Params.eta = 3.8;              % Forcings of equilibrium CO2 doubling (GAMS fco22x)
Params.M_AT_Base = 588;        % Base atm carbon concentration (GAMS in FORC(t) eqn)
Params.delta = 0.1;            % Capital depreciation (5 year) (GAMS dk)
Params.gamma = 0.3;            % Capital elasticity in production function (GAMS gama)
Params.theta2 = 2.8;           % Exponent of control cost function (GAMS expcost2)
Params.alpha = 1.45;           % Elasticity of marginal utility of consumption (GAMS elasmu)
Params.rho = 0.015;            % Initial rate of social time preference per year (GAMS prstp)
Params.xi1 = 0.098;            % Climate equation coefficient for upper level (GAMS c1)
Params.xi2 = 5/3.666;          % Conversion factor from GtC to CtCO2

% Climate Model Diffusion Parameters
%   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
%   DICE2013R-mc uses standard row-column indexing.
phi11 = 0.862962206896552;
phi12 = 0.008624;
phi21 = 0.025;
phi22 = 0.975;

Params.Phi_T = [phi11 phi12; phi21 phi22];

% Carbon Cycle Model Diffusion Parameters 
%   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
%   DICE2013R-mc uses standard row-column indexing.
zeta11 = 0.912;
zeta21 = 0.088;
zeta12 = (588/1350)*zeta21;
zeta22 = 1 - zeta12 - 0.00250;
zeta32 = 0.00250;
zeta23 = zeta32*(1350/10000);
zeta33 = 1 - zeta23;

Params.Phi_M = [zeta11 zeta12 0; zeta21 zeta22 zeta23; 0 zeta32 zeta33];

