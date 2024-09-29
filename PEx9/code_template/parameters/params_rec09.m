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

function params = params_rec09()
%% Control Parameters
params.ctrl.name = 'robust performance learning-based MPC';
params.ctrl.N = 8;
params.ctrl.Q = 1*eye(2);
params.ctrl.R = 100;

%% Estimated System Parameters
% system paramters
params.sys.linear = true; %[bool]
params.sys.d = 4; %[1/s^2]
params.sys.c = 1.2; %[1/s]
params.sys.l = 1.3; %[m]
params.sys.g = 9.81; %[m/s^2]
params.sys.T_s = 0.1; %[s]
params.sys.alpha_s = 10*pi/180; % [rad], linearization point for linear system
params.sys.p = 1; % Number of estimated parameters
params.sys.theta = 0; %[1/s]
params.sys.n=2;
params.sys.m=1;

% system matrices A = A_delta{1} + d*A_delta{2}
params.sys.A_delta{1}=[1, params.sys.T_s;
    -params.sys.T_s*params.sys.d+ params.sys.T_s*params.sys.g/params.sys.l*cos(params.sys.alpha_s), 1-params.sys.T_s*params.sys.c];
params.sys.A_delta{2}=[0,0;-params.sys.T_s, 0];
params.sys.B_delta{1} = [0;params.sys.T_s];
params.sys.C = eye(params.sys.n);
params.sys.D = zeros(params.sys.n);
% state constraints
params.sys.A_x = [1,0; -1,0; 0,1; 0,-1];
params.sys.b_x = [30*pi/180; 30*pi/180; 30*pi/180; 30*pi/180]; %[rad]
% input constraints
params.sys.A_u = [1;-1];
params.sys.b_u = [5;5]; %[1/s^2]
% noise description
params.sys.A_w = [1,0;-1,0;0,1;0,-1];
params.sys.b_w = 0.01*ones(4,1);
params.sys.noiseMean = zeros(params.sys.n,1);
params.sys.noiseCovariance = 0.1*eye(params.sys.n);
params.sys.noiseArgs = {params.sys.noiseMean, params.sys.noiseCovariance, params.sys.A_w, params.sys.b_w};
params.sys.generateNoise = @gaussian_trunc;
% parameter description
params.sys.A_theta = [];
params.sys.b_theta = [];
params.sys.paramVariance = 0;

% uncertainty description
params.sys.d_min = 3;
params.sys.d_max = 5;
params.sys.betaMean = zeros(1,1);
params.sys.betaVar = 10;

%% True System Parameters
% system paramters
params.sys_true.linear = true; %[bool]
params.sys_true.d = params.sys.d_min; %[1/s^2]
params.sys_true.k_unknown = false;
params.sys_true.c = 1.2; %[1/s]
params.sys_true.c_unknown = false;
params.sys_true.l = 1.3; %[m]
params.sys_true.g = 9.81; %[m/s^2]
params.sys_true.T_s = 0.1; %[s]
params.sys_true.alpha_s = 10*pi/180; % [rad], linearization point for linear system
params.sys_true.p = 0; % Number of estimated parameters
params.sys_true.n=2;
params.sys_true.m=1;

% system matrices A = A_delta{1} + d*A_delta{2}
params.sys_true.A=[1, params.sys_true.T_s;
    -params.sys_true.T_s*params.sys_true.d+ params.sys_true.T_s*params.sys_true.g/params.sys_true.l*cos(params.sys_true.alpha_s), 1-params.sys_true.T_s*params.sys_true.c];
params.sys_true.B = [0;params.sys_true.T_s];
params.sys_true.C = eye(params.sys_true.n);
params.sys_true.D = zeros(params.sys_true.n);

% state constraints
params.sys_true.A_x = [1,0; -1,0; 0,1; 0,-1];
params.sys_true.b_x = [30*pi/180; 30*pi/180; 30*pi/180; 30*pi/180]; %[rad]
% input constraints
params.sys_true.A_u = [1;-1];
params.sys_true.b_u = [5;5]; %[1/s^2]
% noise description
params.sys_true.A_w = [1,0;-1,0;0,1;0,-1];
params.sys_true.b_w = 0.01*ones(4,1);
params.sys_true.noiseMean = zeros(params.sys.n,1);
params.sys_true.noiseCovariance = 0.1*eye(params.sys.n);
params.sys_true.noiseArgs = {params.sys.noiseMean, params.sys.noiseCovariance, params.sys.A_w, params.sys.b_w};
params.sys_true.generateNoise = @gaussian_trunc;
% parameter description
params.sys_true.A_theta = [];
params.sys_true.b_theta = [];
params.sys_true.paramVariance = 0;

%% Simulation Parameters
params.sim.nrSteps = 200;
params.sim.nrTraj = 1;
params.sim.x_0 = [0*pi/180;0*pi/180]; %[rad]

%% Plot Parameters
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 1;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end
