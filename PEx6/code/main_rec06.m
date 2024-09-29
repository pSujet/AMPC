%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2022 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reset everything and configure path
clear all;
close all;
clc;
% fix random number generator
rng(6);

%% Robust Constraint Tightening MPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get parameters and define system and controller
params = get_params('params_rec06');
sys = LinearSystem(params.sys);
ctrl = constraint_tightening_RMPC(sys, params.ctrl);

%% Exercise 1a
% Implement compute_robust_tightening in the constraint_tightening_RMPC
% class. Then, run the code below.

% compute infinite horizon LQR controller
[K, P] = dlqr(sys.A, sys.B, params.ctrl.Q, params.ctrl.R);
K = -K;

% compute disturbance reachable sets and corresponding tightenings
F = ctrl.compute_robust_tightening(K); % TODO: Implement the tightening for RMPC

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
figIdx = 1:11;
plot_rec06(figIdx(1), params.sim.nrSteps, X_tight, U_tight, ...
           sys.X, sys.U, F, params.plot)

%% Exercise 1b
% Implement the robust constraint tightening MPC problem in constraint_tighening_RMPC.
% Then, we simulate the closed-loop system for nrSteps time steps and nrTraj
% trajectories starting in x_0.

sys = LinearSystem(params.sys);
ctrl = constraint_tightening_RMPC(sys, params.ctrl, K, P);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0;

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
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

% plot states and inputs
plot_x('state-time', figIdx(2), x, sys.X, params.plot);
plot_x('state-state', figIdx(3), x, sys.X, params.plot);
plot_u(figIdx(4), u, sys.U, params.plot);


%% Stochastic Constraint Tightening MPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get parameters and define system and controller
p = 0.8; % probability level for stochastic backoff term
sys = LinearSystem(params.sys);
ctrl = constraint_tightening_SMPC(sys, params.ctrl,K,P,p);

%% Exercise 2a
% Implement compute_robust_tightening and compute_stochastic_tightening in
% the constraint_tightening_SMPC class. Then, run the code below.

% compute infinite horizon LQR controller
[K, P] = dlqr(sys.A, sys.B, params.ctrl.Q, params.ctrl.R);
K = -K;

% compute disturbance reachable sets and corresponding tightenings
F = ctrl.compute_robust_tightening(K);  % TODO: Implement the tightening for SMPC
% compute stochastic backoff term
[Fw_x, Fw_u] = ctrl.compute_stochastic_tightening(p);  % TODO: Implement the tightening for SMPC

% compute tightened state constraints
X_tight{1} = sys.X;
X_tight{2} = sys.X - Fw_x;
for i=3:length(F)
    X_tight{i} = sys.X - (sys.A + sys.B*K)*F{i-1} - Fw_x;
end

% compute tightened input constraints
U_tight{1} = sys.U;
U_tight{2} = sys.U - Fw_u;
for i=3:length(F)
    U_tight{i} = sys.U - K*(sys.A + sys.B*K)*F{i-1} - Fw_u;
end

% compute DRS used in stochastic constraint tightening MPC
F_s{1} = Polyhedron();
F_s{2} = Polyhedron();
for i=3:length(F)
    F_s{i} = (sys.A + sys.B*K)*F{i-1};
end

% visualize tightenings
plot_rec06(figIdx(5), params.sim.nrSteps, X_tight, U_tight, ...
           sys.X, sys.U, F_s, params.plot)
       
%% Exercise 2b
% Implement the stochastic constraint tightening MPC problem in constraint_tighening_SMPC.
% Then, we simulate the closed-loop system for nrSteps time steps and nrTraj
% trajectories starting in x_0.

sys = LinearSystem(params.sys);
ctrl = constraint_tightening_SMPC(sys, params.ctrl,K,P,p);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0;

% allocate state and input trajectories
x = zeros(nrSteps+1,size(x_0,1),nrTraj); 
u = zeros(nrSteps,nrTraj);
x(1,:,:) = repmat(x_0,[1,1,nrTraj]);

% control-loop
for i=1:nrTraj
    for j=1:nrSteps
        [sol, ~, errmsg] = ctrl.solve(x(j,:,i)', {}, 0);
        % we need this additional logic in the stochastic case, since the
        % MPC problem is allowed to become infeasible and we need a
        % recovery mechanism.
       if errmsg
           error(errmsg);
       end
       u(j,i) = sol(1);
       x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

% plot states and inputs
plot_x('state-time', figIdx(6), x, sys.X, params.plot);
plot_x('state-state', figIdx(7), x, sys.X, params.plot);
plot_u(figIdx(8), u, sys.U, params.plot);

%% Comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 3a: Compare region of attractions of RMPC and SMPC
% Run the code below for different choices of p and observe how the two
% region of attractions change. Please close Figure 9 before rerunning the
% code below, this avoids plotting over the previous result.

p = 1; % probability level for stochastic backoff term
sys = LinearSystem(params.sys);

% compute infinite horizon LQR controller
[K, P] = dlqr(sys.A, sys.B, params.ctrl.Q, params.ctrl.R);
K = -K;

% define robust and stochastic constraint tightening MPC controllers
r_ctrl = constraint_tightening_RMPC(sys, params.ctrl,K,P);
s_ctrl = constraint_tightening_SMPC(sys, params.ctrl,K,P,p);

% grid state constraints
x_0_grid = sys.X.grid(20)';

% allocate region of attractions (RoA)
r_RoA = nan(size(x_0_grid,2),1);
s_RoA = nan(size(x_0_grid,2),1);

% check grid for feasibility
for i=1:size(x_0_grid,2)
    % solve RMPC
    try
        [r_sol, ~, ~] = r_ctrl.solve(x_0_grid(:,i), {}, 0);
        r_RoA(i) = logical(r_sol); 
    catch
        r_RoA(i) = 0;
    end
    
    % solve SMPC
    try
        [s_sol, ~, ~] = s_ctrl.solve(x_0_grid(:,i), {}, 0);
        s_RoA(i) = logical(s_sol);
    catch
        s_RoA(i) = 0;
    end
end

% compute RoA of RMPC
r_RoA = Polyhedron(x_0_grid(:,logical(r_RoA))');
r_RoA.computeHRep();
% compute RoA of SMPC
s_RoA = Polyhedron(x_0_grid(:,logical(s_RoA))');
s_RoA.computeHRep();

% plot RoAs
figure(figIdx(9));
s_RoA.plot('color', 'k', 'linestyle', '--', 'linewidth', 0.2, 'alpha', 0.2);
hold on;
r_RoA.plot('color', 'k', 'linestyle', '--', 'linewidth', 0.2, 'alpha', 0.6);
hold on;
sys.X.plot('wire', true, 'linestyle', '-', 'linewidth', 1.5);

legend("RoA SMPC", "RoA RMPC", "State Constraints");

%% Let's change some parameters
params = get_params('params_rec06');
params.sim.nrSteps = 5;
params.sim.nrTraj = 150;
params.sim.x_0 = [25*pi/180;-22*pi/180];

%% Exercise 3b: Compare closed-loop trajectories of RMPC and SMPC
% Run the code below for different choices of p and observe how the
% closed-loop trajectories and constraint violations change. Please close
% Figures 10 & 11 before rerunning the code below, this avoids plotting
% over the previous result.

p = 0.3; % probability level for stochastic backoff term

% 1. compute closed-loop trajectories of RMPC
% basic definitions
sys = LinearSystem(params.sys);
ctrl = constraint_tightening_RMPC(sys, params.ctrl,K,P);
nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0;

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
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

% create phase plot
plot_x('state-state', figIdx(10), x, sys.X, params.plot);

% 2. compute closed-loop trajectories of SMPC
% basic definitions
sys = LinearSystem(params.sys);
ctrl = constraint_tightening_SMPC(sys, params.ctrl,K,P,p);

nrSteps = params.sim.nrSteps;
nrTraj = params.sim.nrTraj;
x_0 = params.sim.x_0;

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
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
        %end
    end
end

% create phase plot
plot_x('state-state', figIdx(11), x, sys.X, params.plot);
