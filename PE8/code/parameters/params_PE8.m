%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022, Alexandre Didier and Jérôme Sieber, ETH Zurich,
% {adidier,jsieber}@ethz.ch
%
% All rights reserved.
%
% This code is only made available for students taking the advanced MPC 
% class in the fall semester of 2022 (151-0371-00L) and is NOT to be 
% distributed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = params_PE8()
%% Control Parameters
params.ctrl.name = 'Model Predictive Safety Filter';
params.ctrl.N = 30;

%% System Parameters
% system dimensions
n = 2; m = 1;
params.sys.n = n;
params.sys.m = m;

% system paramters
params.sys.T_s = 0.1; %[s]
k = 8; g = 9.81; l = 1.3; c = 1;
params.sys.A = [1, params.sys.T_s;
                params.sys.T_s*(-k+g/l), 1 - params.sys.T_s*c];
params.sys.B = [0;params.sys.T_s];
params.sys.C = eye(n);
params.sys.D = zeros(n,m);

% state constraints
params.sys.A_x = [1,0; -1,0; 0,1; 0,-1];
params.sys.b_x = [45*pi/180; 30*pi/180; 30*pi/180; 30*pi/180]; %[rad]
% input constraints
params.sys.A_u = [1;-1];
params.sys.b_u = [2;2]; %[1/s^2]
% noise description
params.sys.A_w = [];
params.sys.b_w = [];

params.sys.generateNoise = @zero;
params.sys.noiseArgs = {n};

%% Simulation Parameters
params.sim.nrSteps = 100;
params.sim.nrTraj = 1;
params.sim.x_0 = [10*pi/180; -20*pi/180]; %[rad]

%% Plot Parameters
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 1;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end
