function plot_gams_verification(u, x, SCC, J, x0, Params)
%
% Plot the solution of DICE2013R_mc and
% the reference solution obtained with GAMS
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
gamma       = Params.gamma;
theta2      = Params.theta2;

% Translate time indices into years
end_year = 2010 + (N+1)*5;
years = 2010:5:end_year;

% -------------------------------------------------------------------------
%   Load GAMS results from comma separated value file
gams_results = csvread('Gams_Results.csv',0,1,[0 1 14 60]);
g_years = 2010:5:end_year-10;

% -------------------------------------------------------------------------
%   Plot Exogenous Signals

figure(5), title('Exogenous Signals')

subplot(3,2,1) 
stem(years,L), title('Population'), xlabel('years'), hold on
plot(g_years,gams_results(6,:),'r','LineWidth',2)
legend('DICE2013R\_mc','DICE2013R - GAMS Outputs')
temp = axis; axis([2010 2305 6000 13000]);

subplot(3,2,2) 
stem(years,A_TFP), title('Total Factor Productivity'), xlabel('years'), hold on
plot(g_years,gams_results(2,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,3) 
stem(years,sigma), title('Carbon Intensity'), xlabel('years'), hold on
plot(g_years,gams_results(3,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,4) 
stem(years(1:end-1),theta1), title('\theta_1'), xlabel('years'), hold on
plot(g_years,gams_results(4,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,5) 
stem(years(1:end-1),E_Land), title('Land Use Change Emissions'), xlabel('years'), hold on
plot(g_years,gams_results(1,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,6) 
stem(years(1:end-1),F_EX), title('Exogenous Forcings'), xlabel('years'), hold on
plot(g_years,gams_results(5,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

% -------------------------------------------------------------------------
% Plot Endogenous States

figure(6), title('Endogenous States')

subplot(3,2,1), 
stem(years(1:end-1),x(1,:)), ylabel('TAT'), xlabel('years'), hold on
plot(g_years,gams_results(7,:),'r','LineWidth',2)
legend('DICE2013R\_mc','DICE2013R - GAMS Outputs')
temp = axis; axis([2010 2305 temp(3) temp(4)+0.5]);

subplot(3,2,2), 
stem(years(1:end-1),x(2,:)), ylabel('TLO'), xlabel('years'), hold on
plot(g_years,gams_results(8,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,3) 
stem(years(1:end-1),x(3,:)), ylabel('MAT'), xlabel('years'), hold on
plot(g_years,gams_results(9,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,4) 
stem(years(1:end-1),x(4,:)), ylabel('MUP'), xlabel('years'), hold on
plot(g_years,gams_results(10,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,5) 
stem(years(1:end-1),x(5,:)), ylabel('MLO'), xlabel('years'), hold on
plot(g_years,gams_results(11,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

subplot(3,2,6) 
stem(years(1:end-1),x(6,:)), ylabel('K'), xlabel('years'), hold on
plot(g_years,gams_results(12,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);
 
% -------------------------------------------------------------------------
% Plot Control Values

figure(7)

subplot(2,1,1) 
stem(years(1:end-2),u(1,:)), title('Emissions Control Rate'), xlabel('years'), hold on
plot(g_years,gams_results(13,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);
legend('DICE2013R\_mc','DICE2013R - GAMS Outputs')

subplot(2,1,2) 
stem(years(1:end-2),u(2,:)), title('Savings Rate'), xlabel('years'), hold on
plot(g_years,gams_results(14,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 temp(3) temp(4)]);

% -------------------------------------------------------------------------
% Plot Social Cost of Carbon

figure(8)

stem(years(1:end-2),SCC), hold on
xlabel('years')
ylabel('2005 USD')
title('Social Cost of Carbon')
plot(g_years,gams_results(15,:),'r','LineWidth',2)
temp = axis; axis([2010 2305 0 400]);
legend('DICE2013R\_mc','DICE2013R - GAMS Outputs')

% -------------------------------------------------------------------------
% Calculate Optimal Welfare from GAMS results
% Compute optimal trajectories from GAMS optimal control
xgams = zeros(7,N+1);
Gross_Economic_Outputg = zeros(1,N);
Emissionsg = zeros(1,N);
Damagesg = zeros(1,N);
Net_Economic_Outputg = zeros(1,N);
Per_Cap_Consumptiong = zeros(1,N);
xgams(:,1) = x0(1:7);
ugams = gams_results(13:14,:);
for t = 1:N
    xgams(:,t+1) = dice_dynamics(xgams(:,t),ugams(:,t),t,Params);
    Gross_Economic_Outputg(t) = A_TFP(t) * (xgams(6,t)^gamma) * ((L(t)/1000)^(1-gamma));
    Emissionsg(t) = sigma(t) * (1-ugams(1,t)) * Gross_Economic_Outputg(t) + E_Land(t);
    Damagesg(t) = 1 - 0.00267*(xgams(1,t)^2);
    Net_Economic_Outputg(t) = (Damagesg(t) - theta1(t)*(ugams(1,t)^theta2)) * Gross_Economic_Outputg(t);
    Per_Cap_Consumptiong(t) = (1000*(1-ugams(2,t))/L(t)) * Net_Economic_Outputg(t);
end

display('============================================================')
format long; display(['J (Matlab):   ', num2str(J)])      
format long; display(['J (GAMS)  :   ', num2str(-xgams(7,N+1))])    
display('============================================================')
