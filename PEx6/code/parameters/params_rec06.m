%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022, ETH Zurich, {adidier, jsieber}@ethz.ch
%
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2022 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = params_rec06()
%% Control Parameters
params.ctrl.name = 'constraint tightening SMPC';
params.ctrl.N = 6;
params.ctrl.Q = eye(2);
params.ctrl.R = 1000;

%% System Parameters
% system dimensions
n = 2; m = 1;
params.sys.n = n;
params.sys.m = m;
params.sys.dt = 0.15;

% linear system
k = 0; g = 9.81; l = 1.3; c = 0.5;
params.sys.A = [1, params.sys.dt;
                params.sys.dt*(-k+g/l), 1 - params.sys.dt*c];
params.sys.B = [0;params.sys.dt];
params.sys.C = eye(n);
params.sys.D = zeros(n,m);

% state constraints
params.sys.A_x = [1,0; -1,0; 0,1; 0,-1];
params.sys.b_x = [25*pi/180; 25*pi/180; 25*pi/180; 25*pi/180]; %[rad]
% input constraints
params.sys.A_u = [1;-1];
params.sys.b_u = [4;4]; %[1/s^2]
% noise description
params.sys.A_w = [1,0; -1,0; 0,1; 0,-1];
params.sys.b_w = [2*pi/180; 2*pi/180; 2.2*pi/180; 2.2*pi/180]; %[rad, rad, rad/s, rad/s]
% noise distribution
params.sys.generateNoise = @uniform;
params.sys.noiseArgs = {params.sys.A_w, params.sys.b_w};

%% Simulation Parameters
params.sim.nrSteps = 30;
params.sim.nrTraj = 50;
params.sim.x_0 = [20*pi/180;-5*pi/180]; %[rad]

%% Plot Parameters
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 0.2;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end