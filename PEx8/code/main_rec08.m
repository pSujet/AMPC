%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2022 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reset everything and configure path
setup;
% fix random number generator
rng(8);
% get parameters
params = get_params('params_rec08');

%% Initialize data
% Data set
x_USA=4.4;
y_USA=11.7;
x_UK=7.6;
y_UK=19.4;
x_SW=6.6;
y_SW=30.0;
x_AT=8.1;
y_AT=25.1;
x_FR=4.3;
y_FR=10.7;
x_JP=1.2;
y_JP=2.2;

%% Problem Set Ex 1
% Compute the posterior distribution p(theta | x_USA, y_USA)
% TODO: Complete the code for the constructor in estimators/BLR.m


%%% TODO: Initialise the data and the linear features
%%% Change the following lines to compute the posterior given the data of
%%% the USA
%%% Hint: features should be a function handle [@(x)] returning x given x
data_init.x = x_USA;
data_init.y = y_USA;

% linear features
features = @(x) x;

% Prior for linear features
theta_p = 5;
Sigma_p = 4;

linear_BLR = BLR(params.sys, theta_p, Sigma_p, data_init, features);

theta_USA=linear_BLR.theta;
Sigma_USA=linear_BLR.Sigma;

disp('Posterior distribution given the data point for the USA:')
disp(['Mean: ', num2str(theta_USA)])
disp(['Variance: ', num2str(Sigma_USA)])

% plot
figIdx = 1:4;
plot_rec08(figIdx(1), data_init, linear_BLR, [], params.plot)

%% Problem Set Ex 2
% Compute the recursive BLR update given x_UK and y_UK
% TODO: Complete the update_estimate function in estimators/BLR.m

%%% TODO: Initialise the data and the linear features 
%%% Change the following lines to update the BLR given the new data of the
%%% UK
data_update.x=x_UK;
data_update.y=y_UK;

linear_BLR.update_estimate(data_update);

theta_post = linear_BLR.theta;
Sigma_post = linear_BLR.Sigma;

disp('Posterior distribution given the data point for the USA and the UK:')
disp(['Mean: ', num2str(theta_post)])
disp(['Variance: ', num2str(Sigma_post)])

% plot
plot_rec08(figIdx(2), struct('x', [data_init.x, data_update.x], 'y', [data_init.y, data_update.y]), linear_BLR, [], params.plot)

%% Problem Set Ex 3
% Predict the amount of nobel laureates for Switzerland based on the
% provided data
% TODO: Complete the predict metod in estimators/BLR.m

x_CH=8.8;
[y_CH, var_CH] = linear_BLR.predict(x_CH);

disp('Prediction of the number of swiss Nobel laureates per 10 million population using BLR:')
disp(['Value: ', num2str(y_CH)])
disp(['Variance: ', num2str(var_CH)])

% plot
plot_rec08(figIdx(3), struct('x', [data_init.x, data_update.x], 'y', [data_init.y, data_update.y]), linear_BLR, [x_CH; y_CH], params.plot)

%% Problem Set Ex 4
% Compute the BLR with nonlinear features with zero mean and identity prior
data.x=[x_USA, x_UK, x_SW, x_AT, x_FR, x_JP];
data.y=[y_USA, y_UK, y_SW, y_AT, y_FR, y_JP];

%%% TODO: Initialise the nonlinear features
%%% Change the following lines for the BLR with nonlinear features
%%% Features should be a function handle returning the stacked vector x and 
%%% x_i^2.

% Features
n=2;
features = @(x)[x; x.^2];

% Prior
theta_p = zeros(n,1);
Sigma_p = 0.5*eye(n);

nonlinear_BLR = BLR(params.sys, theta_p, Sigma_p, data, features);

% Prediction
x_CH=8.8;
[y_CH, var_CH] = nonlinear_BLR.predict(x_CH);

disp('Prediction of the number of swiss Nobel laureates per 10 million population using BLR with nonlinear features:')
disp(['Prediction: ', num2str(y_CH)]);
disp(['Variance: ', num2str(var_CH)]);

% plot
plot_rec08(figIdx(4), data, nonlinear_BLR, [x_CH; y_CH], params.plot)