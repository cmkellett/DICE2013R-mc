function plot_results(u, x, SCC, J, Params)
%
% Plot the solution of DICE2013R_mc
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

N           = Params.N;
sigma       = Params.sigma;
A_TFP       = Params.A_TFP;
L           = Params.L;
theta1      = Params.theta1;
F_EX        = Params.F_EX;
E_Land      = Params.E_Land;

% Translate time indices into years
end_year = 2010 + (N+1)*5;
years = 2010:5:end_year;

% -------------------------------------------------------------------------
%   Plot Exogenous Signals

figure(1), title('Exogenous Signals')

subplot(3,2,1) 
stem(years,L), title('Population'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,2) 
stem(years,A_TFP), title('Total Factor Productivity'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,3) 
stem(years,sigma), title('Carbon Intensity'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,4) 
stem(years(1:end-1),theta1), title('\theta_1'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,5) 
stem(years(1:end-1),E_Land), title('Land Use Change Emissions'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,6) 
stem(years(1:end-1),F_EX), title('Exogenous Forcings'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

% -------------------------------------------------------------------------
% Plot Endogenous States

figure(2), title('Endogenous States')

subplot(3,2,1), 
stem(years(1:end-1),x(1,:)), ylabel('TAT'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,2), 
stem(years(1:end-1),x(2,:)), ylabel('TLO'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,3) 
stem(years(1:end-1),x(3,:)), ylabel('MAT'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,4) 
stem(years(1:end-1),x(4,:)), ylabel('MUP'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,5) 
stem(years(1:end-1),x(5,:)), ylabel('MLO'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,6) 
stem(years(1:end-1),x(6,:)), ylabel('K'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);
 
% -------------------------------------------------------------------------
% Plot Control Values

figure(3)

subplot(2,1,1) 
stem(years(1:end-2),u(1,:)), title('Emissions Control Rate'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(2,1,2) 
stem(years(1:end-2),u(2,:)), title('Savings Rate'), xlabel('years'), hold on
temp = axis; axis([2010 2305 temp(3) temp(4)]);

% -------------------------------------------------------------------------
% Plot Social Cost of Carbon

figure(4)

stem(years(1:end-2),SCC)
xlabel('years')
ylabel('2005 USD')
title('Social Cost of Carbon')
temp = axis; axis([2010 2305 temp(3) temp(4)]);

display('============================================================')
format long; display(['Optimal Welfare:   ', num2str(J)])      
display('============================================================')
