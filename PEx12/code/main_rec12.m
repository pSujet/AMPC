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

%% reset everything and configure path
clear all;
close all;
clc;

% fix random number generator
rng(12);

figIdx = 1:10;

% get parameters and initialize system
params = get_params('params_rec12');
sys = NonlinearSystem(params.sys);

%% Define simple 1D example for unknown function h
h = @(x) -x.^2 .* sin(5 .* pi .* x).^6; % domain x \in [0,1]
h_derivative = @(x) -2.*x.*sin(5.*pi.*x).^6 + 30*pi.*x.^2.*(sin(5*pi.*x).^5).*cos(5.*pi.*x);
% true minimum
x_min = 0.9015; y_min = -0.8113;

% plot and derivative
figure(figIdx(1));
title('Simple 1D Example for h');
X = 0:0.001:1;
xlabel('x');
yyaxis left
plot(X, h(X)); hold on; ylim([-0.9, 0.9]); ylabel('y'); hold on;
yyaxis right
plot(X, h_derivative(X)); ylabel('dy/dx'); hold on; plot([0,1], [0,0], 'k'); hold on; ylim([-25, 25]);

%% Sample a few values of unknown function to fit GP
nrSamples = 10;
data.x = rand(1,nrSamples);
data.y = h(data.x);

% define GP
% low noise covariance since we don't have noise
params.gp.noiseCovariance = sqrt(0.01);
% squared exponential kernel
lambda = 1; sigma = 0.05;
params.gp.kernel = @(x,y) lambda^2*exp(-abs(x'-y).^2/(2*sigma^2));

% initialize BO
bo = BO(data, params.gp);

% plot fitted GP
plot_rec12(figIdx(2), bo, h, [], [], params.plot);

%% TODO: Implement acquisition function in BO class
% sample acquisition function
nrSamples = 100;
X = 0:1/nrSamples:1;
[theta, y_hat] = bo.acquisition_fn(X);

% plot updated GP and acquisition function
plot_rec12(figIdx(3), bo, h, y_hat, X, params.plot);

%% Sample suggested parameter and add new data to GP
data.x = theta;
data.y = h(theta);

% add new data
bo.add_data(data);

% sample acquisition function
X = 0:1/nrSamples:1;
[theta, y_hat] = bo.acquisition_fn(X);

% plot updated GP and acquisition function
plot_rec12(figIdx(3), bo, h, y_hat, X, params.plot);

%% Run BO for a few iterations
nrItr = 10;

for i=1:nrItr
    data.x = theta;
    data.y = h(theta);

    bo.add_data(data);

    X = 0:1/nrSamples:1;
    theta = bo.acquisition_fn(X);
end

% plot updated GP
plot_rec12(figIdx(4), bo, h, [], [], params.plot);

%% Get best estimate from fitted GP
nrSamples = 200;
X = 0:1/nrSamples:1;
[theta, y, conf] = bo.get_estimate(X);

disp(['Best estimated parameter: ', num2str(theta), '. True parameter: ',  num2str(x_min)]);
disp(['Predicted value [95% confidence interval]: ', num2str(y), ' [', num2str(conf),']. True value at minimum: ',  num2str(y_min)]);

%% Run PID on nonlinear segway system
%  trajectories starting in x_0
ctrl = PID(sys, params.ctrl);

% simulate system with defined controller
[x, u] = simulate(sys, ctrl, params.sim);

% plot trajectories
plot_x('state-time', figIdx(5), x, sys.X, params.plot);
plot_x('state-state', figIdx(6), x, sys.X, params.plot);
plot_u(figIdx(7), u, sys.U, params.plot);

%% Run PID on nonlinear segway system again with different PID gains
% set PID gains
params.ctrl.Kp = 8;
params.ctrl.Ki = 6;
params.ctrl.Kd = 0.9;
ctrl = PID(sys, params.ctrl);

% simulate system with defined controller
[x, u] = simulate(sys, ctrl, params.sim);

% plot trajectories
plot_x('state-time', figIdx(8), x, sys.X, params.plot);
plot_x('state-state', figIdx(9), x, sys.X, params.plot);
plot_u(figIdx(10), u, sys.U, params.plot);