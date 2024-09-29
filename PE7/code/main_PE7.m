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

%% reset everything, configure path, and define necessary variables
clear all;
close all;
clc;

% fix random number generator
rng(12);

figIdx = 1:11;

% get parameters and initialize system
params = get_params('params_PE7');
sys = NonlinearSystem(params.sys);

% Define simple 1D example for unknown function h
% objective function
h = @(x) -x.^2 .* sin(5 .* pi .* x).^6; % domain x \in [0,1]
h_derivative = @(x) -2.*x.*sin(5.*pi.*x).^6 + 30*pi.*x.^2.*(sin(5*pi.*x).^5).*cos(5.*pi.*x);
% true minimum
x_min = 0.9015; y_min = -0.8113;

%% PART 1: BO on simple 1D example
% Sample a few values of unknown function to fit GP
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
plot_PE7_2D(figIdx(1), bo, h, [], [], params.plot);

%% TODO: Implement acquisition function (LCB) in BO class
% sample acquisition function
nrSamples = 250;
Q = 0:1/nrSamples:1;
[q, q_hat] = bo.acquisition_fn(Q); % TODO: Implement LCB

plot_PE7_2D(figIdx(2), bo, h, q_hat, Q, params.plot);

%% Sample suggested parameter and add new data to GP
% Use this cell to visualize how BO explores the parameter space by
% rerunning this cell for a couple of times. Make sure to close the figure
% after each run to avoid overlay.
% Please note that the lower confidence bound (LCB) might not align with
% the shaded area depending on which value you chose for beta!

data.x = q;
data.y = h(q);

% add new data
bo.add_data(data);

% sample acquisition function
Q = 0:1/nrSamples:1;
[q, q_hat] = bo.acquisition_fn(Q);

% plot updated GP and acquisition function
plot_PE7_2D(figIdx(2), bo, h, q_hat, Q, params.plot);

%% Run BO for a few iterations
nrItr = 30;

for i=1:nrItr
    data.x = q;
    data.y = h(q);
    
    % add data to GP in BO instance
    bo.add_data(data);
    
    % sample acquisition function
    Q = 0:1/nrSamples:1;
    q = bo.acquisition_fn(Q);
end

% plot updated GP
plot_PE7_2D(figIdx(3), bo, h, [], [], params.plot);

%% Get best estimate from fitted GP
nrSamples = 200;
Q = 0:1/nrSamples:1;
[q, y, conf] = bo.get_estimate(Q);

disp(['Best estimated parameter: ', num2str(q), '. True parameter: ',  num2str(x_min)]);
disp(['Predicted value [95% confidence interval]: ', num2str(y), ' [', num2str(conf),']. True value at minimum: ',  num2str(y_min)]);

%% PART 2: BO for PID tuning
% we only optimize the proportional and integral gains to allow easy
% visualization of the fitted GP, therefore the range of Kd is set to a
% singleton

% define ranges for PID gains Kp, Ki, Kd
params.ctrl.range = [0, 10;     % Kp
                     0, 10;     % Ki
                     0.9, 0.9]; % Kd

%% Run PID on nonlinear segway system once
% Let's see how well some random PID gains perform
ctrl = PID(sys, params.ctrl);

% simulate system with defined controller
[x, u] = simulate(sys, ctrl, params.sim);

% plot trajectories
plot_x('state-time', figIdx(4), x, sys.X, params.plot);
plot_x('state-state', figIdx(5), x, sys.X, params.plot);
plot_u(figIdx(6), u, sys.U, params.plot);

%% Randomly sample a few data points to fit GP
% generate a few parameter vectors randomly
nrSamples = 50;
Q = params.ctrl.range(:,1) + (params.ctrl.range(:,2) - params.ctrl.range(:,1)).*rand(3,nrSamples);

% allocate output data vector
y = zeros(nrSamples,1); 

for i=1:nrSamples
    % set PID gains
    params.ctrl.Kp = Q(1,i);
    params.ctrl.Ki = Q(2,i);
    params.ctrl.Kd = Q(3,i);
    ctrl = PID(sys, params.ctrl);

    % simulate system with chosen PID gains
    [x, u] = simulate(sys, ctrl, params.sim);
    
    %%% TODO: Implement a cost function that incentivises PID gains that
    %%% can perform a swing-up and keep the pendulum (aka segway) in an
    %%% upright position. Additionally, the cost should penalize jittering
    %%% in the control action (please also see the exercise sheet).
    %%% Hint: In this exercise you have a lot of freedom, many choices of
    %%% cost functions will work. However, as a starting point have a look
    %%% at the recitation (e.g. the recording) to get some ideas.
    % --- start inserting code here ---
    y(i) = sum(x(end-50:end,1).^2) + sum(x(end-50:end,2).^2) + sum(u.^2);
    % --- stop inserting code here ---
end

%% Fit GP by instanciating BO
data.x = Q;
data.y = y';

%%% TODO: Implement the squared exponential kernel for the GP used in the
%%% BO class.
% --- start inserting code here ---
% low noise covariance since we don't have noise
params.gp.noiseCovariance = sqrt(0.01);
% squared exponential kernel
lambda = 2; sigma = 0.7; % you can change these hyperparameters if you want, however the provided ones should work fine
params.gp.kernel = @(x,y) lambda^2*exp(-(reshape(vecnorm(reshape(reshape(x,[],1)-repmat(y,size(x,2),1),3,[])),[],size(reshape(x,[],1)-repmat(y,size(x,2),1),2))).^2/(2*sigma^2)); % implement kernel here
% --- stop inserting code here ---

% initialize GP
bo = BO(data, params.gp);

%% visualize fitted GP
% the resulting plot is a 3D plot, however MATLAB might only show it in 2D
% initially. You can rotate the plot in any direction using the rotating
% tool.
nrSamples = 1000;
Q = params.ctrl.range(:,1) + (params.ctrl.range(:,2) - params.ctrl.range(:,1)).*rand(3,nrSamples);

% plot fitted GP
plot_PE7_3D(figIdx(7), bo, Q, params);

%% Apply BO for some iterations

nrItr = 200;
for i=1:nrItr
    % sample acquisition function
    nrSamples = 100;
    Q = params.ctrl.range(:,1) + ( params.ctrl.range(:,2) - params.ctrl.range(:,1)).*rand(3,nrSamples);
    q = bo.acquisition_fn(Q);

    % set PID gains
    params.ctrl.Kp = q(1);
    params.ctrl.Ki = q(2);
    params.ctrl.Kd = q(3);
    ctrl = PID(sys, params.ctrl);
    
    % simulate system with chosen PID gains
    [x, u] = simulate(sys, ctrl, params.sim);
    
    %%% TODO:
    % Copy the cost function from above here and correctly set up the data
    % struct to pass to the BO class
    % --- start inserting code here ---
    data.x = [q(1);q(2);q(3)];   
    data.y = sum(x(end-50:end,1).^2) + sum(x(end-50:end,2).^2) + sum(u.^2);
    % --- stop inserting code here ---

    % add new data point to GP in BO instance
    bo.add_data(data);
end

%% visualize the resulting GP and visualize resulting PID controller
% the resulting GP plot is a 3D plot, however MATLAB might only show it in
% 2D initially. You can rotate the plot in any direction using the rotating
% tool.

nrSamples = 1000;
Q = params.ctrl.range(:,1) + (params.ctrl.range(:,2) - params.ctrl.range(:,1)).*rand(3,nrSamples);
% get best estimate for PID gains
[q, y, conf] = bo.get_estimate(Q);

% set PID gains
params.ctrl.Kp = q(1);
params.ctrl.Ki = q(2);
params.ctrl.Kd = q(3);
ctrl = PID(sys, params.ctrl);

% simulate system with chosen PID gains
[x, u] = simulate(sys, ctrl, params.sim);

% plot trajectories
plot_x('state-time', figIdx(8), x, sys.X, params.plot);
plot_x('state-state', figIdx(9), x, sys.X, params.plot);
plot_u(figIdx(10), u, sys.U, params.plot);
% plot fitted GP
plot_PE7_3D(figIdx(11), bo, Q, params);

disp('Best estimated PID gains:');
disp(['K_P: ', num2str(q(1)), ', K_I: ',  num2str(q(2)), ', K_D: ',  num2str(q(3))]);
disp(['Predicted cost [95% confidence interval]: ', num2str(y), ' [', num2str(conf),'].']);