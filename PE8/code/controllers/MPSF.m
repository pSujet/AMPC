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

classdef MPSF < Controller
    %MPSF Safety Filter Class
    %   Class for safety filter 
    
    properties
        
    end
    
    methods
        function obj = MPSF(sys, params, P, solver)
            %MPSF Construct an instance of this class
            %   Construct safety filter Class and initialize solver
            %%% Parse inputs %%%
            switch nargin
                case 3
                    solver = '';
                    
                case 4
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            obj@Controller(sys, params);
            
            %%% Exercise 5: Initialize the MPC Controller           
            
            % SDP Variables
            x = sdpvar(sys.n, obj.params.N+1);
            x_0 = sdpvar(sys.n,1);
            u = sdpvar(sys.m, obj.params.N);
            u_L = sdpvar(sys.m, 1);
            
            %%% TODO: Implement the model predictive safety filter
            %%% optimisation problem. The input matrix P corresponds to the
            %%% ellipsoid {x| x'*P*x <= 1} which you computed in Exercise
            %%% 2.
            %%% Using sedumi will lead to numerical issues in this
            %%% exercise, so we use fmincon.
            %%% If you plan on using Mosek, make sure you pass 'mosek-socp'
            %%% in the solver argument instead of just 'mosek'.
            % --- start inserting code here ---
            % objectives
            objective = norm(u-u_L)^2;            
            % constraints
            constraints = [x(:,1) == x_0];
            for i = 1:obj.params.N
                constraints = [constraints, x(:,i+1) == sys.A*x(:,i)+sys.B*u(:,i)];
                constraints = [constraints, sys.X.A*x(:,i)<=sys.X.b];
                constraints = [constraints, sys.U.A*u(:,i)<=sys.U.b];
            end
            constraints = [constraints, x(:,obj.params.N+1)'*P*x(:,obj.params.N+1)<=1];            
            % --- stop inserting code here ---
            
            % Setup solver object
            obj.prob=optimizer(constraints, objective, sdpsettings('solver',solver), {x_0, u_L}, {u(:,1)});
            
        end
    end
end

