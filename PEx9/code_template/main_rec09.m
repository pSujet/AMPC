cd c:\gurobi1000\win64\matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% Reset everything and configure path
setup
% fix random number generator
rng(9);
% get parameters for this recitation
params = get_params('params_rec09');

%% Exercise 2: Extend exercise of previous recitation to GPs 
% Compute the posterior distribution p(theta | x_CH, y_CH)
% Complete the implementation of the kernel handles below.

% initialize given data points
x_USA=4.4;
y_USA=11.7;
x_UK=7.6;
y_UK=19.4;
x_SW=6.6;
y_SW=30.0;
data.x =[x_USA, x_UK, x_SW];
data.y =[y_USA, y_UK, y_SW];

% define noise distribution
params.noiseCovariance = sqrt(2);


%%% TODO %%%
% Implement the quadratic kernel.
%  Play around with a different kernel, e.g., the squared exponential
% kernel. Also change the hyperparameters and observe how the prediction
% changes.
% --- start inserting here ---            
% quadratic kernel
kernel_quad = @(x,y)x.^2'*y.^2;

% squared exponential kernel
lambda = 10; sigma = 0.1;
kernel_se = @(x,y) lambda*exp(-0.5*abs(x'-y).^2/(sigma^2));

% initialize GP
gp = GP(params, data, kernel_se); % change the kernel function here

% --- stop inserting here ---
%%%


% predict
x_CH=8.8;
[y_CH, var_CH] = gp.predict(x_CH);

disp('Prediction of the number of Swiss Nobel laureates per 10 million population using a GP:')
disp(['Value: ', num2str(y_CH)])
disp(['Variance: ', num2str(var_CH)])

% Plot the PredictiveDistribution
figIdx = 1:6;

figure(figIdx(1)); hold on;
title('Gaussian Process Predictive Distribution Mean')
x_plot=(4:0.1:10); [mean_plot, var_plot] = gp.predict(x_plot);
plot(x_plot, mean_plot+3*sqrt(diag(var_plot)), 'k--', 'LineWidth',1.5);
plot(x_plot, mean_plot-3*sqrt(diag(var_plot)), 'k--', 'LineWidth',1.5);
p1 = fill([x_plot'; flip(x_plot')], [mean_plot+3*sqrt(diag(var_plot)); flip(mean_plot-3*sqrt(diag(var_plot)))], [220,220,220]/256, 'FaceAlpha', 0.8, 'EdgeAlpha',0);
p2 = plot(x_plot,mean_plot,'k','LineWidth',1.5);
p3 = plot(data.x,data.y,'bx','MarkerSize',10);
p4 = plot(x_CH,y_CH,'r*','MarkerSize',15);
xlabel('x'); ylabel('y'); grid();
legend([p3,p2,p1(1),p4],{'Data','Mean','Uncertainty','Prediction'},'Location','nw');
set(gcf,'position',[100,100,0.7*params.plot.width,1.4*params.plot.height],'color','white')

%% Exercise 3: Robust Performance Learning-based MPC
% Compute steady state input

% look up parameters
alpha_s = params.sys.alpha_s;
d = params.sys.d;
g = params.sys.g;
l = params.sys.l;

% compute steady state input
u_s = alpha_s*(d - g/l*cos(alpha_s));

%% Compute disturbance set

% look up parameters
T_s = params.sys.T_s;
d_hat = params.sys.d;
d_max = params.sys.d_max;
d_min = params.sys.d_min;
B = params.sys.B_delta{1};
A_x = params.sys.A_x;
b_x = params.sys.b_x;

%%% TODO %%%
% Compute the set of all values of g(Delta_x, Delta_u) by using all
% vertices of A_x*Delta_x<=b_x and all vertices of Delta_d
% Hint: You only need to consider the vertices for the second state
% --- start inserting here ---            

g_vertices = [];
% --- stop inserting here ---
%%%

% Disturbance set
g_max = max(g_vertices);
g_min = min(g_vertices);
params.sys.b_w = params.sys.b_w+[0; 0; g_max; -g_min]; % Mikowski Sum

%% visualize computed disturbance set
sys = LinearAffineSystem(params.sys);
sys_true = LinearSystem(params.sys_true);
ctrl = RPLBMPC(sys, params.ctrl);

% compute infinite horizon LQR controller
[K, P] = dlqr(sys.A, sys.B, params.ctrl.Q, params.ctrl.R);
K = -K;

% compute disturbance reachable sets and corresponding tightenings
F = ctrl.compute_robust_tightening(K);

% compute tightened state constraints
X_tight{1} = sys.X;
for i=2:length(F)
    X_tight{i} = sys.X - F{i};
end

% compute tightened input constraints
U_tight{1} = sys.U;
for i=2:length(F)
    U_tight{i} = sys.U - K*F{i};
end

% plot tightenings
plot_rec9(figIdx(2), params.sim.nrSteps, X_tight, U_tight, ...
           sys.X, sys.U, F, params.plot)

%% Simulate Robust Performance Learning-based MPC
% TODO: Implement the RPLMPC problem in the constructor of the RPLBMPC
% class.
u_s_true = params.sys_true.alpha_s*(params.sys_true.d - params.sys_true.g/params.sys_true.l*cos(params.sys_true.alpha_s));
ctrl = RPLBMPC(sys, params.ctrl);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0 - [alpha_s; 0]; % need this to convert initial state to new coordinate frame

% initialize BLR estimator
noise_params.noiseMean = 0;
noise_params.noiseCovariance = 0.1; 
noise_params.noiseArgs = {noise_params.noiseMean, noise_params.noiseCovariance};

betaMean = params.sys.betaMean;
betaVar = params.sys.betaVar;
features = @(x) [x];
data.x =[]; data.y = [];

linear_BLR = BLR(noise_params, betaMean, betaVar, data , features);

% allocate state and input trajectories
x = zeros(nrSteps+1,size(x_0,1),nrTraj); 
u = zeros(nrSteps,nrTraj);
x(1,:,:) = repmat(x_0,[1,1,nrTraj]);
% allocate Delta_d trajectories
Delta_d = zeros(nrSteps+1,nrTraj);

% control-loop
for i=1:nrTraj
    for j=1:nrSteps
        [sol, ~, errmsg] = ctrl.solve(x(j,:,i)', {Delta_d(j,i)}, 0);
        if errmsg
            error(errmsg);
        end
        u(j,i) = sol(1);
        x(j+1,:,i) = sys_true.step(x(j,:,i)', u(j,i)-(u_s_true-u_s));
        
        % parameter estimation
        data.x = T_s*x(j,2,i)-T_s*alpha_s;
        y = x(j+1,:,i)' - sys.A*x(j,:,i)' - sys.B*(u(j,i));
        data.y = y(2);
        linear_BLR.update_estimate(data);
        Delta_d(j+1,i) = linear_BLR.theta;        
    end
end

% plot results
params.plot.color = [0.7216, 0.1490, 0.0039];
p1 = plot_x('state-time', figIdx(3), x + repmat([alpha_s, 0],[nrSteps+1,1,nrTraj]), sys.X, params.plot);
p2 = plot_x('state-state', figIdx(4), x + repmat([alpha_s, 0],[nrSteps+1,1,nrTraj]), sys.X, params.plot);
p3 = plot_u(figIdx(5), u - u_s, sys_true.U, params.plot);

% plot parameter learning
params.plot.color = [0.7216, 0.1490, 0.0039];
figure(figIdx(6)); 
plot(Delta_d,'color',[params.plot.color,params.plot.alpha],'linewidth',params.plot.lw);
hold on; plot([0, nrSteps+2], [-(params.sys.d-params.sys_true.d),-(params.sys.d-params.sys_true.d)], 'k--', 'linewidth', 1.);
xlabel('time'); ylabel('Delta\_d'); xlim([0,nrSteps+2]); ylim([min(min(Delta_d),-(params.sys.d-params.sys_true.d))-0.01, max(max(Delta_d),-(params.sys.d-params.sys_true.d))+0.01]); grid();
hold on; title('Learn Delta\_d using BLR');
set(gcf,'position',[100,100,params.plot.width,2*params.plot.height],'color','white')

%% compare to robust constraint tightening MPC
sys = LinearAffineSystem(params.sys);
ctrl = constraint_tightening_RMPC(sys, params.ctrl, K, P);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0 - [alpha_s; 0]; % need this to convert initial state to new coordinate frame

% allocate state and input trajectories
x = zeros(nrSteps+1,size(x_0,1),nrTraj); 
u = zeros(nrSteps,nrTraj);
x(1,:,:) = repmat(x_0,[1,1,nrTraj]);

% control-loop
for i=1:nrTraj
    for j=1:nrSteps
        [sol, ~, errmsg] = ctrl.solve(x(j,:,i)', {}, 0);
        if errmsg
            error(errmsg);
        end
        u(j,i) = sol(1);
        x(j+1,:,i) = sys_true.step(x(j,:,i)', u(j,i) - (u_s_true-u_s));
    end
end

% plot results/comparison
params.plot.color = [0.392, 0.584, 0.929];
p4 = plot_x('state-time', figIdx(3), x + repmat([alpha_s, 0],[nrSteps+1,1,nrTraj]), sys.X, params.plot); hold on; 
subplot(2,1,1); p_s = plot([1,nrSteps+1], [alpha_s*180/pi, alpha_s*180/pi], 'k--', 'linewidth', 1.2); hold on;
subplot(2,1,2); plot([1,nrSteps+1], [0, 0], 'k--', 'linewidth', 1.2); 
legend([p1(1); p4(1); p_s], {'Robust Performance LBMPC', 'Robust Constraint Tightening MPC', 'Steady State'}, 'Location', 'southeast')
p5 = plot_x('state-state', figIdx(4), x + repmat([alpha_s, 0],[nrSteps+1,1,nrTraj]), sys.X, params.plot);
legend([p2(1); p5(1)], {'Robust Performance LBMPC', 'Robust Constraint Tightening MPC'}, 'Location', 'northeast')
p6 = plot_u(figIdx(5), u - u_s, sys_true.U, params.plot);
legend([p3(1); p6(1)], {'Robust Performance LBMPC', 'Robust Constraint Tightening MPC'}, 'Location', 'northeast')
