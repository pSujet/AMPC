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
setup
% fix random number generator
rng(13);

% get parameters and define system and safety filter
params = get_params('params_PE8');
sys = LinearSystem(params.sys);

% "learning"-based input 
u_L=sin(0.5*pi*(0:params.sys.T_s:params.sys.T_s*params.sim.nrSteps));

figIdx = 1:10;

%% Exercise 1: simulate the open-loop system for nrSteps time steps and nrTraj
% TODO: look at the figures and verify that the constraints are violated 

%  trajectories starting in x_0
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
        u(j,i) = u_L(j);
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

% plot results for "learn
plot_x('state-time', figIdx(1), x, sys.X, params.plot);
plot_x('state-state', figIdx(2), x, sys.X, params.plot);
plot_u(figIdx(3), u, sys.U, params.plot);

%% Exercise 2,3,4: Compute ellipsoidal safe set
%TODO: Implement the invariant set computation in IBSF.m
%      Write the solve function in IBSF.m
sf = IBSF(sys, params.ctrl);

plot_PE8(figIdx(4), params.sim.nrSteps, sys.X, sys.U, sf.P, 1, params.plot)

%% simulate the closed-loop system for nrSteps time steps and nrTraj
%  trajectories starting in x_0
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
        u(j,i) = sf.solve(x(j,:,i)',u_L(j));
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

%% plot results
plot_x('state-time', figIdx(5), x, sys.X, params.plot);
plot_x('state-state', figIdx(6), x, sys.X, params.plot);
plot_PE8(figIdx(6), params.sim.nrSteps, sys.X, sys.U, sf.P, 1, params.plot)
plot_u(figIdx(7), u, sys.U, params.plot);

%% Exercise 5: MPSF
% TODO: Write the MPSF Problem in MPSF.m. Use the ellipsoid computed in
% exercise 2 as the terminal set
% Note: We use fmincon here as sedumi has numerical issues
% Note: If you plan on using Mosek, use the SOCP solver, i.e. 'mosek-socp'
% and not just 'mosek'. This can lead to Yalmip using the LP/QP solver for
% an SOCP if the constraint is defined quadratically!
sf_MPSF = MPSF(sys, params.ctrl, sf.P,'fmincon');

%% simulate the closed-loop system for nrSteps time steps and nrTraj
%  trajectories starting in x_0
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
        [sol, ~, errmsg] = sf_MPSF.solve(x(j,:,i)',{u_L(j)});
        if errmsg
            error(errmsg);
        end
        u(j,i) = sol;
        x(j+1,:,i) = sys.step(x(j,:,i)', u(j,i));
    end
end

%% plot results
plot_x('state-time', figIdx(8), x, sys.X, params.plot);
plot_x('state-state', figIdx(9), x, sys.X, params.plot);
plot_PE8(figIdx(9), params.sim.nrSteps, sys.X, sys.U, sf.P, 1, params.plot)
plot_u(figIdx(10), u, sys.U, params.plot);