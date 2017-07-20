function [E_Industrial, Net_Output, Per_Cap_Consumption, Damages_Fraction, ...
    Atm_Carbon_ppm, Marg_Cost_Abatement] = ...
    compute_auxiliary_quantities(x_opt,Params,nx);
% Computes additional outputs based on the solution of the DICE
%   optimal control problem.
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
% -------------------------------------------------------------------------
%                           Unpack States and Controls
% -------------------------------------------------------------------------

% DICE inputs: mitigation rate and savings rate
mu = x_opt(nx+1,:);
s = x_opt(nx+2,:);

% DICE endogenous states
TATM = x_opt(1,:);
TLO = x_opt(2,:);
MATM = x_opt(3,:);
MUP = x_opt(4,:);
MLO = x_opt(5,:);
K = x_opt(6,:);

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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                           Compute quantities of interest
% -------------------------------------------------------------------------
for i=1:N
    Gross_Economic_Output(i) = A_TFP(i)*(K(i)^gamma)*((L(i)/1000)^(1-gamma));
    Damages(i) = 1 - 0.00267*(TATM(i)^2); 
    
    Damages_Fraction(i) = 1 - Damages(i);
    Net_Output(i) = (Damages(i) - theta1(i)*(mu(i)^theta2))*Gross_Economic_Output(i);
    E_Industrial(i) = sigma(i)*(1-mu(i))*Gross_Economic_Output(i);
    Per_Cap_Consumption(i) = 1000*(1-s(i)) * Net_Output(i) / L(i);
    Atm_Carbon_ppm(i) = MATM(i)/2.13;
    Marg_Cost_Abatement(i) = 344 * (0.975^(i-1)) * (mu(i)^1.8);
end