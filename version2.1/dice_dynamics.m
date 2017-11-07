function [f, varargout] = dice_dynamics(x,u,i,Params,varargin)
%
% function [f, varargout] = dice_dynamics(x,u,i,Params,varargin)
%
% Evolution of endogenous state variables, including emissions and 
%   consumption for compution of social cost of carbon.
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

% number of additional input and output arguments
noutextra = max(nargout)-1;
ninextra  = max(nargin) -4; 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                               States
% -------------------------------------------------------------------------
T_AT = x(1); 
T_LO = x(2);
T = [T_AT; T_LO];

M_AT = x(3);
M_UP = x(4);
M_LO = x(5);

M = [M_AT; M_UP; M_LO];

K = x(6);

J = x(7);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                               Controls
% -------------------------------------------------------------------------

% Mitigation Rate
mu = u(1);

% Savings Rate
s = u(2);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                           Unpack Parameters
% -------------------------------------------------------------------------

N           = Params.N;
sigma       = Params.sigma;
A_TFP       = Params.A_TFP;
L           = Params.L;
theta1      = Params.theta1;
F_EX        = Params.F_EX;
E_Land      = Params.E_Land;
eta         = Params.eta;
M_AT_Base   = Params.M_AT_Base;
delta       = Params.delta;
gamma       = Params.gamma;
theta2      = Params.theta2;
alpha       = Params.alpha;
rho         = Params.rho;
xi1         = Params.xi1;
xi2         = Params.xi2;
Phi_T       = Params.Phi_T;
Phi_M       = Params.Phi_M;
zeta11      = Phi_M(1,1);
zeta12      = Phi_M(1,2);

a2          = Params.a2;
a3          = Params.a3;

% Named functional quantities
Gross_Economic_Output = A_TFP(i)*(K^gamma)*((L(i)/1000)^(1-gamma));
Damages = 1 - a2*(T_AT^a3); 
Net_Economic_Output = (Damages - theta1(i)*(mu^theta2))*Gross_Economic_Output;

Emissions_rhs = sigma(i)*(1-mu)*Gross_Economic_Output + E_Land(i);
Consumption_rhs = (1-s) * Net_Economic_Output; 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                     Handling Variables of Interest
% -------------------------------------------------------------------------
switch ninextra
    case 1 % use symbolic variable v to determine SCC
        v =varargin{1};
        Emissions   = v(1); % Emission
        Consumption = v(2); % Consumption
    case 0 % standard computation
        Emissions = Emissions_rhs;
        Consumption =  Consumption_rhs;
end

Radiative_Forcing = xi1*(eta*log((zeta11*M_AT + zeta12*M_UP + xi2*Emissions)/M_AT_Base)/log(2) + F_EX(i+1));

% ---------------- Climate ------------------------
T_NEXT = Phi_T * T + [Radiative_Forcing; 0];

% ---------------- Carbon Cycle -------------------
M_NEXT = Phi_M * M + [(5/3.666)*Emissions; 0; 0];

% ---------------- Economy ------------------------
K_NEXT = (1 - delta)^5 * K + 5 * Net_Economic_Output * s;

% ---------------- Objective ----------------------
J_NEXT = J - L(i)*(((1000/L(i)*Consumption)^(1-alpha) - 1)/((1-alpha))-1)/(1+rho)^(5*(i-1));

f = [T_NEXT; M_NEXT; K_NEXT; J_NEXT];

% Extra output arguments needed for SCC computation
if noutextra == 1
   varargout{1} = [Emissions_rhs; Consumption_rhs]; 
end

