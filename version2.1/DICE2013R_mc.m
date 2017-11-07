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
%
%   This code provides a bare bones replication of the GAMS DICE2013R
%   (http://aida.wss.yale.edu/~nordhaus/homepage/Web-DICE-2013-April.htm)
%   for use in MATLAB with CasADi.  See the accompanying documentation
%   for more details.
%
%
% -------------------------------------------------------------------------
% Release notes
%
%  - 1.0 initial release requires casadi-matlabR2014b-v3.1.1 or v.3.0.0 
%  - 1.1 fixes compatability with casadi-matlabR2014b-v3.2.1
%  - 2.0 adds parameters for DICE2016R-091916ap, with switch to choose
%       between 2013R (110513_vanilla) and 2016R parameter values
%  - 2.1 fixes minor parameter errors in 2016R and time indexing errors
%       in sigma and the radiative forcing
% -------------------------------------------------------------------------

clc
clear all
close all

% import casadi name space
import casadi.*

%% ========================================================================
% define data of the DICE Optimal Control Problem (OCP)
% =========================================================================
Param_set = 2013; % Choose DICE 2013, 2016, or 1 for custom values 
Params = set_DICE_parameters(Param_set); % get DICE parameters
N = Params.N;

switch Param_set
    case 1
        x0 = [0.8 0.0068 830.4 1527 10010 135]'; % initial condition 
    % IT IS RECOMMENDED THAT CHANGES BE LIMITED TO THE ABOVE!!!
    % Initial conditions below correspond to DICE2013R and DICE2016R.
    case 2013
        x0 = [0.8 0.0068 830.4 1527 10010 135]'; % initial condition 
    case 2016
        x0 = [0.85 0.0068 851 460 1740 223]'; % initial condition
end
x0 = [x0;0]; % extra state for objective computation

nx = length(x0); % # of states
nu = 2; % # of inputs
nv = 2; % # of variables of interest (consumption & emission) needed for efficient computation of SCC 

%% ========================================================================
% construct guess for decision variables of Nonlinear Program (NLP)
% =========================================================================
u_guess(1,:) = [0 0.5*ones(1,N-3) 0 0 0];
u_guess(2,:) = [0.5*ones(1,N) 0];

x_guess = zeros(nx, N+1);

for i = 1:N
    if i == 1
        x_guess(:,1) = x0(1:nx);
    end
  [fi, hi] = dice_dynamics(x_guess(:,i),u_guess(:,i),i,Params);
  x_guess(:, i+1) = fi;
  v_guess(:, i)   = hi;
end
v_guess(:, N+1)   = hi;

xuv_guess = [x_guess(:);u_guess(:);v_guess(:)];

%% ========================================================================
% construct bounds on decision variables of NLP
% =========================================================================
u_LB = [zeros(1,N+1); [zeros(1,N-10) Params.optlsrv*ones(1,10) 0]];
u_UB = [[Params.miu0 1.0*ones(1,28) 1.2*ones(1,N-28)]; ones(1,N+1)];
u_LB = u_LB(:,1:N+1);
u_UB = u_UB(:,1:N+1);

x_LB = zeros(nx, N+1);
x_LB(end,:) = -inf; % no lower bound on objective 
x_UB = inf*ones(nx, N+1); % no upper bound on states 

v_LB = -inf*ones(nv, N+1);
v_UB =  inf*ones(nv, N+1);

xuv_LB = [x_LB(:); u_LB(:); v_LB(:)];
xuv_UB = [x_UB(:); u_UB(:); v_UB(:)];

%% ========================================================================
% define NLP
% =========================================================================
tic
% allocate CASADI variables
x      = SX.sym('x', nx, N+1); % states + value 
u      = SX.sym('u', nu, N+1); % inputs
v      = SX.sym('v', nv, N+1); % variables of interest, i.e. emissions & consumption
eq_con = SX.sym('eq_con', nx+nv, N+1); % (nx+nv) * N+1 constraints for the dynamics and variables of interest

% loop over dynamics 
for  i = 1:N
    if  i == 1
       eq_con(1:nx,1) = x(1:nx,1) - x0(1:nx); % x(:,1) = x0;    
    end
    % equality constraints for dynamics
    [fi, hi] = dice_dynamics(x(1:nx,i),u(:,i),i,Params,v(:,i));
    eq_con(1:nx,i+1) = x(1:nx,i+1) - fi; 
    
    % equality constraints assigning emmissions and consumption 
    eq_con(nx+1:end, i) = (v(:,i) - hi);
    if i == N
     eq_con(nx+1:end, i+1) = v(:,i+1);    % dummy constraint at i = N+1
    end
    
end

% define the objective (Mayer term)
obj = ((5 * Params.scale1 * x(end, N+1)) - Params.scale2);

% define nlp variables
nlp = struct('x', [x(:);u(:);v(:)], 'f', obj, 'g', eq_con(:));

% options for IPOPT
opts = struct;
opts.ipopt.max_iter    = 3000;
opts.ipopt.print_level = 3;%0,3
opts.print_time        = 0;
opts.ipopt.acceptable_tol =1e-12;
opts.ipopt.acceptable_obj_change_tol = 1e-12;
opts.ipopt.mu_strategy                 = 'adaptive';

% create IPOPT solver object
solver = nlpsol('solver', 'ipopt', nlp);
t_CON = toc;

% solve the NLP
tic

lbg = zeros(nx+nv, N+1);
lbg(end-1:end, end) =  0;% % bound for dummy constraint
ubg = zeros(nx+nv, N+1); 
ubg(end-1:end, end) =  0; % bound for dummy constraint

res = solver( 'x0' , xuv_guess,...      % solution guess
              'lbx', xuv_LB,...         % lower bound on x
              'ubx', xuv_UB,...         % upper bound on x
              'lbg', 0,...              % lower bound on g
              'ubg', 0);                % upper bound on g
t_NLP = toc;

% print computation times to screen
display('============================================================')
display(['Time to construct NLP:     ', num2str(t_CON)])
display(['Time to solve NLP:         ', num2str(t_NLP)])
display('============================================================')

%% ========================================================================
% assign and plot results
% =========================================================================
x_res = full(res.x);
% reshape solution to matrix variables
v_opt = x_res(end-length(v(:))+1:end);
v_opt = reshape(v_opt, nu, N+1);

u_opt = x_res(end-length(v(:))-length(u(:))+1:end-length(v(:)));
u_opt = reshape(u_opt, nu, N+1);
u_opt = u_opt(:,1:end-1); % last value is not meaningful

x_opt = x_res(1:end-length(v(:))-length(u(:)));
x_opt = reshape(x_opt, nx, N+1);

% Compute auxiliary outputs
[E_Industrial, Net_Output, Per_Cap_Consumption, Damages_Fraction, ...
    Atm_Carbon_ppm, Marg_Cost_Abatement] = ...
    compute_auxiliary_quantities(x_opt,u_opt,Params,nx);

% Optimal Welfare
J = -x_opt(nx,end);

% DICE inputs: mitigation rate and savings rate
mu = u_opt(1,:);
s  = u_opt(2,:);

% DICE endogenous states
TATM = x_opt(1,:);
TLO = x_opt(2,:);
MATM = x_opt(3,:);
MUP = x_opt(4,:);
MLO = x_opt(5,:);
K = x_opt(6,:);
x_opt = x_opt(1:nx-1,:);

% Marginals and Social Cost of Carbon
lam  = full(res.lam_g);
lam  = reshape(lam, nx+nv, N+1);
lamE = lam(end-1,1:end);
lamC = lam(end,1:1:end);
if Param_set == 2016
    SCC = -1000*[lamE(1:end-1)]./[.00001+lamC(1:end-1)];
else
    SCC = -1000*[lamE(1:end-1)]./[lamC(1:end-1)];
end


% Tidy Workspace.  Comment out the following two commands to maintain all
%   variables in workspace.
clear_list = {'eq_con','fi','hi','i','lam','lbg','nlp','nu','nv','nx', ...
    'obj','opts','res','solver','t_CON','t_NLP','u','u_guess','u_LB', ...
    'u_UB','ubg','v','v_guess','v_LB','v_UB','x_guess','x_LB','x_UB', ...
    'xuv_guess','xuv_LB','xuv_UB'};
clear(clear_list{:});

% Default plot results
plot_results(u_opt,x_opt,SCC,J,Params,Param_set);

% Uncomment the following line to generate plots comparing the output
%   of DICE2013R_mc and the outputs obtained via the original GAMS
%   code for DICE2013R.  Note, this requires the file GAMS_Results.csv.
% plot_gams_verification(u_opt,x_opt,SCC,J,x0,Params,Param_set);
