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

function params = params_PE7()
%% Control Parameters
n = 2; m = 1;
params.ctrl.name = 'PID';
params.ctrl.setpoint = zeros(n,1);
params.ctrl.ss_input = zeros(m,1);
params.ctrl.Kp = 8;
params.ctrl.Ki = 2;
params.ctrl.Kd = 0.05;

%% System Parameters
% system dimensions
params.sys.n = n;
params.sys.m = m;
params.sys.Ts = 0.05;

% nonlinear system dynamics
k = 3; c = 0; g = 9.81; l = 1.3;
params.sys.f = @(x,u) [x(1) + params.sys.Ts*x(2);
                       x(2) + params.sys.Ts*(-k*x(1) - c*x(2) + sin(x(1))*g/l + u)];
params.sys.h = @(x,u) x;

% state constraints
params.sys.A_x = [];%[1,0; -1,0; 0,1; 0,-1];
params.sys.b_x = [];%[pi/4; pi/4; pi/3; pi/3]; %[rad]
% input constraints
params.sys.A_u = [1;-1];
params.sys.b_u = [2;2]; %[1/s^2]
% noise description
params.sys.A_w = [];
params.sys.b_w = [];
% noise distribution
params.sys.generateNoise = @zero;
params.sys.noiseArgs = {n};

%% Simulation Parameters
params.sim.nrSteps = 200;
params.sim.nrTraj = 1;
params.sim.x_0 = [180*pi/180;0]; %[rad]

%% Plot Parameters
params.plot.show = true;
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 1;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end