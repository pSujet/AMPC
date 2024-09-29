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

function [x,u] = simulate(sys, ctrl, params)
%SAMPLE Sample a system given a controller and return
%       state and input trajectories

    %%% Parse input arguments %%%
    switch nargin
        case 3
            
        otherwise
            error('Wrong number of inputs!')
    end
    %%%%%%%%%%%%%%%%%%%
    
    %  trajectory starting in x_0
    nrSteps = params.nrSteps;
    x_0 = params.x_0;

    % allocate state and input trajectories
    x = zeros(nrSteps+1,size(x_0,1)); 
    u = zeros(nrSteps,1);
    x(1,:) = repmat(x_0,[1,1]);

    % control-loop
    for i=1:nrSteps
        [sol, ~] = ctrl.solve(x(i,:)');
        
        u(i) = sol;
        x(i+1,:) = sys.step(x(i,:)', u(i));
    end
end

