function Params = set_DICE_parameters(Param_set)
%
% Params = set_DICE_parameters(Param_set)
%
% Param_set = release version (2013, 2016) or 1 for custom values
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

switch Param_set
    case 1
    % -------------------------------------------------------------------------
    % Constants
    Params.N = 60;                 % Horizon length
    Params.BaseYear = 2010;        % Initial year
    Params.eta = 3.8;              % Forcings of equilibrium CO2 doubling (GAMS fco22x)
    Params.M_AT_Base = 588;        % Base atm carbon concentration (GAMS in FORC(t) eqn)
    Params.delta = 0.1;            % Capital depreciation (5 year) (GAMS dk)
    Params.gamma = 0.3;            % Capital elasticity in production function (GAMS gama)
    Params.theta2 = 2.8;           % Exponent of control cost function (GAMS expcost2)
    Params.a2 = 0.00267;           % Damage multiplier
    Params.a3 = 2;                 % Damage exponent
    Params.alpha = 1.45;           % Elasticity of marginal utility of consumption (GAMS elasmu)
    Params.rho = 0.015;            % Initial rate of social time preference per year (GAMS prstp)
    Params.xi1 = 0.098;            % Climate equation coefficient for upper level (GAMS c1)
    Params.xi2 = 5/3.666;          % Conversion factor from GtC to CtCO2
    Params.scale1 = 0.016408662;   % Utility multiplier in cost function
    Params.scale2 = 3855.106895;   % Utility offset in cost function
    
    % Climate Model Diffusion Parameters
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    c1 = Params.xi1;
    c3 = 0.088;
    c4 = 0.025;
    fco22x = Params.eta;
    t2xco2 = 2.9;

    phi11 = 1-c1*((fco22x/t2xco2) + c3);
    phi12 = c1*c3;
    phi21 = c4;
    phi22 = 1-c4;

    Params.Phi_T = [phi11 phi12; phi21 phi22];

    % Carbon Cycle Model Diffusion Parameters 
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    b12 = 0.088;
    b23 = 0.0025;
    mateq = 588;
    mueq = 1350;
    mleq = 10000;

    zeta11 = 1 - b12;
    zeta21 = b12;
    zeta12 = (mateq/mueq)*zeta21;
    zeta22 = 1 - zeta12 - b23;
    zeta32 = b23;
    zeta23 = zeta32*(mueq/mleq);
    zeta33 = 1 - zeta23;

    Params.Phi_M = [zeta11 zeta12 0; zeta21 zeta22 zeta23; 0 zeta32 zeta33];

    % -------------------------------------------------------------------------
    % Exogenous Signal Constants
    Params.L0 = 6838;           % Initial population
    Params.La = 10500;          % Asymptotic population
    Params.lg = 0.134;          % Population growth rate

    Params.EL0 = 3.3;           % Initial land use emissions
    Params.deltaEL = 0.2;       % Land use emissions decrease rate

    Params.A0 = 3.8;            % Initial Total Factor Productivity (TFP)
    Params.ga = 0.079;          % Initial TFP rate
    Params.deltaA = 0.006;      % TFP increase rate

    Params.pb = 344;            % Initial backstop price
    Params.deltaPB = 0.025;     % Decline rate of backstop price

    e0 = 33.61;                 % Initial emissions
    q0 = 63.69;                 % Initial global output
    Params.miu0 = 0.039;        % Initial mitigation rate
    Params.sigma0 = e0/(q0*(1-Params.miu0)); % Calculated initial emissions intensity
    Params.gsigma = 0.01;       % Emissions intensity base rate
    Params.deltasigma = 0.001;  % Decline rate of emissions intensity

    Params.f0 = 0.25;           % Initial forcings of non-CO2 GHGs
    Params.f1 = 0.7;            % Forcings of non-CO2 GHGs in 2100
    Params.tforce = 18;         % Slope of non-CO2 GHG forcings

    Params.optlsrv = (Params.delta + .004)/(Params.delta ...
        + 0.004*Params.alpha + Params.rho) * Params.gamma;
    
% -------------------------------------------------------------------------
% IT IS STRONGLY RECOMMENDED THAT YOU DO NOT CHANGE PARAMETERS BELOW !!!
% Parameter values for DICE2013R-110513_vanilla and DICE2016R-091916ap 
% are provided below.
    
    case 2013
    % -------------------------------------------------------------------------
    % Constants
    Params.N = 60;                 % Horizon length
    Params.BaseYear = 2010;        % Initial year
    Params.eta = 3.8;              % Forcings of equilibrium CO2 doubling (GAMS fco22x)
    Params.M_AT_Base = 588;        % Base atm carbon concentration (GAMS in FORC(t) eqn)
    Params.delta = 0.1;            % Capital depreciation (5 year) (GAMS dk)
    Params.gamma = 0.3;            % Capital elasticity in production function (GAMS gama)
    Params.theta2 = 2.8;           % Exponent of control cost function (GAMS expcost2)
    Params.a2 = 0.00267;           % Damage multiplier
    Params.a3 = 2;                 % Damage exponent
    Params.alpha = 1.45;           % Elasticity of marginal utility of consumption (GAMS elasmu)
    Params.rho = 0.015;            % Initial rate of social time preference per year (GAMS prstp)
    Params.xi1 = 0.098;            % Climate equation coefficient for upper level (GAMS c1)
    Params.xi2 = 5/3.666;          % Conversion factor from GtC to CtCO2
    Params.scale1 = 0.016408662;   % Utility multiplier in cost function
    Params.scale2 = 3855.106895;   % Utility offset in cost function

    % Climate Model Diffusion Parameters
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    c1 = Params.xi1;
    c3 = 0.088;
    c4 = 0.025;
    fco22x = Params.eta;
    t2xco2 = 2.9;

    phi11 = 1-c1*((fco22x/t2xco2) + c3);
    phi12 = c1*c3;
    phi21 = c4;
    phi22 = 1-c4;

    Params.Phi_T = [phi11 phi12; phi21 phi22];

    % Carbon Cycle Model Diffusion Parameters 
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    b12 = 0.088;
    b23 = 0.0025;
    mateq = 588;
    mueq = 1350;
    mleq = 10000;

    zeta11 = 1 - b12;
    zeta21 = b12;
    zeta12 = (mateq/mueq)*zeta21;
    zeta22 = 1 - zeta12 - b23;
    zeta32 = b23;
    zeta23 = zeta32*(mueq/mleq);
    zeta33 = 1 - zeta23;

    Params.Phi_M = [zeta11 zeta12 0; zeta21 zeta22 zeta23; 0 zeta32 zeta33];

    % -------------------------------------------------------------------------
    % Exogenous Signal Constants
    Params.L0 = 6838;           % Initial population
    Params.La = 10500;          % Asymptotic population
    Params.lg = 0.134;          % Population growth rate

    Params.EL0 = 3.3;           % Initial land use emissions
    Params.deltaEL = 0.2;       % Land use emissions decrease rate

    Params.A0 = 3.8;            % Initial Total Factor Productivity (TFP)
    Params.ga = 0.079;          % Initial TFP rate
    Params.deltaA = 0.006;      % TFP increase rate

    Params.pb = 344;            % Initial backstop price
    Params.deltaPB = 0.025;     % Decline rate of backstop price

    e0 = 33.61;                 % Initial emissions
    q0 = 63.69;                 % Initial global output
    Params.miu0 = 0.039;        % Initial mitigation rate
    Params.sigma0 = e0/(q0*(1-Params.miu0)); % Calculated initial emissions intensity
    Params.gsigma = 0.01;       % Emissions intensity base rate
    Params.deltasigma = 0.001;  % Decline rate of emissions intensity

    Params.f0 = 0.25;           % Initial forcings of non-CO2 GHGs
    Params.f1 = 0.7;            % Forcings of non-CO2 GHGs in 2100
    Params.tforce = 18;         % Slope of non-CO2 GHG forcings

    Params.optlsrv = (Params.delta + .004)/(Params.delta ...
        + 0.004*Params.alpha + Params.rho) * Params.gamma;
    
    case 2016
    % -------------------------------------------------------------------------
    % Constants
    Params.N = 100;                % Horizon length
    Params.BaseYear = 2015;        % Initial year
    Params.eta = 3.6813;           % Forcings of equilibrium CO2 doubling (GAMS fco22x)
    Params.M_AT_Base = 588;        % Base atm carbon concentration (GAMS in FORC(t) eqn)
    Params.delta = 0.1;            % Capital depreciation (5 year) (GAMS dk)
    Params.gamma = 0.3;            % Capital elasticity in production function (GAMS gama)
    Params.theta2 = 2.6;           % Exponent of control cost function (GAMS expcost2)
    Params.a2 = 0.00236;           % Damage multiplier
    Params.a3 = 2;                 % Damage exponent
    Params.alpha = 1.45;           % Elasticity of marginal utility of consumption (GAMS elasmu)
    Params.rho = 0.015;            % Initial rate of social time preference per year (GAMS prstp)
    Params.xi1 = 0.1005;           % Climate equation coefficient for upper level (GAMS c1)
    Params.xi2 = 5/3.666;          % Conversion factor from GtC to GtCO2
    Params.scale1 = 0.030245527;   % Utility multiplier in cost function
    Params.scale2 = 10993.704;     % Utility offset in cost function

    % Climate Model Diffusion Parameters
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    c1 = Params.xi1;
    c3 = 0.088;
    c4 = 0.025;
    fco22x = Params.eta;
    t2xco2 = 3.1;

    phi11 = 1-c1*((fco22x/t2xco2) + c3);
    phi12 = c1*c3;
    phi21 = c4;
    phi22 = 1-c4;

    Params.Phi_T = [phi11 phi12; phi21 phi22];

    % Carbon Cycle Model Diffusion Parameters 
    %   Note Vanilla DICE2013R uses nonstandard column-row indexing, where
    %   DICE2013R-mc uses standard row-column indexing.  The diffusion
    %   parameters are calculated based on other parameters and we include
    %   those here as they change between 2013 and 2016 versions.
    b12 = 0.12;
    b23 = 0.007;
    mateq = 588;
    mueq = 360;
    mleq = 1720;

    zeta11 = 1 - b12;
    zeta21 = b12;
    zeta12 = (mateq/mueq)*zeta21;
    zeta22 = 1 - zeta12 - b23;
    zeta32 = b23;
    zeta23 = zeta32*(mueq/mleq);
    zeta33 = 1 - zeta23;

    Params.Phi_M = [zeta11 zeta12 0; zeta21 zeta22 zeta23; 0 zeta32 zeta33];

    % -------------------------------------------------------------------------
    % Exogenous Signal Constants
    Params.L0 = 7403;           % Initial population
    Params.La = 11500;          % Asymptotic population
    Params.lg = 0.134;          % Population growth rate

    Params.EL0 = 2.6;           % Initial land use emissions
    Params.deltaEL = 0.115;     % Land use emissions decrease rate

    Params.A0 = 5.115;          % Initial Total Factor Productivity (TFP)
    Params.ga = 0.076;          % Initial TFP rate
    Params.deltaA = 0.005;      % TFP increase rate

    Params.pb = 550;            % Initial backstop price
    Params.deltaPB = 0.025;     % Decline rate of backstop price

    e0 = 35.85;                 % Initial emissions
    q0 = 105.5;                 % Initial global output
    Params.miu0 = 0.03;         % Initial mitigation rate
    Params.sigma0 = e0/(q0*(1-Params.miu0)); % Calculated initial emissions intensity
    Params.gsigma = 0.0152;     % Emissions intensity base rate
    Params.deltasigma = 0.001;  % Decline rate of emissions intensity

    Params.f0 = 0.5;             % Initial forcings of non-CO2 GHGs
    Params.f1 = 1.0;             % Forcings of non-CO2 GHGs in 2100
    Params.tforce = 17;          % Slope of non-CO2 GHG forcings
    
    Params.optlsrv = (Params.delta + .004)/(Params.delta ...
        + 0.004*Params.alpha + Params.rho) * Params.gamma;
end
    
% -------------------------------------------------------------------------
% Exogenous Signals
N = Params.N;

Params.sigma = zeros(1,N+1);
Params.L = zeros(1,N+1);
Params.A_TFP = zeros(1,N+1);
Params.E_Land = zeros(1,N+1);
Params.F_EX = zeros(1,N+1);
Params.theta1 = zeros(1,N+1);

sigma(1)    = Params.sigma0;    % Initial sigma (GAMS sig0)
L(1)        = Params.L0;        % Initial population in millions (GAMS pop0)
A_TFP(1)    = Params.A0;        % Initial TFP (GAMS a0)

for i=1:N+1
    sigma(i+1) = sigma(i) * exp(-Params.gsigma * (((1-Params.deltasigma)^5)^(i-1)) * 5);
    L(i+1) = L(i) * (Params.La/L(i))^Params.lg;
    A_TFP(i+1) = A_TFP(i) / (1 - Params.ga * exp(-Params.deltaA * 5 * (i-1)));
    E_Land(i) = Params.EL0*((1-Params.deltaEL)^(i-1));
    F_EX(i) = Params.f0 + min(Params.f1-Params.f0, ...
        ((Params.f1-Params.f0)/Params.tforce)*(i-1));
    theta1(i) = sigma(i) * (Params.pb/(1000*Params.theta2)) * (1-Params.deltaPB)^(i-1);
end

Params.sigma = sigma;
Params.L = L;
Params.A_TFP = A_TFP;
Params.E_Land =E_Land;
Params.F_EX = F_EX;
Params.theta1 = theta1;
