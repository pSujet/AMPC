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

classdef RPLBMPC < Controller
    %RPLBMPC RPLBMPC Controller Class
    %   Construct and solve robust performance LBMPC Problem
    
    properties
        solver % solver name
        
        % since we need to compute the RPI set, we additionally store the
        % system information
        sys % System object
    end
    
    methods
        function obj = RPLBMPC(sys, params, K, P, solver)
            %RPLBMPC Construct an instance of this class
            %   Construct MPC Class and initialize solver
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    [K, P] = dlqr(sys.A, sys.B, params.Q, params.R);
                    K = -K;
                    solver = 'quadprog';
                    
                case 3
                    error('Wrong number of inputs!')
                    
                case 4
                    solver = 'quadprog';
                    
                case 5
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%
            
            % call super class constructor
            obj@Controller(sys, params);
            
            % store system information
            obj.sys = sys;
            
            % initialize solver name
            obj.solver = solver;
            
            %%% Initialize the RPLBMPC Controller
            
            %%% Initialize the Robust Performance LBMPC Controller
            % --- compute robust tightening ---
            F = obj.compute_robust_tightening(K);
            
            % --- compute terminal set ---
            X_f = sys.compute_MRPI(K);
            
            % look up horizon length
            N = obj.params.N;
            
            % look up state and input dimensions
            n = obj.sys.n;
            m = obj.sys.m;
            p = obj.sys.p;
            
            %look up extra parameters
            T_s = obj.sys.params.T_s;
            alpha_s = obj.sys.params.alpha_s;
            
            % look up system matrices
            A = obj.sys.A;
            B = obj.sys.B;
            
            % look up cost matrices
            Q = obj.params.Q;
            R = obj.params.R;
            
            % look up state, input, and disturbance sets
            X = obj.sys.X;
            U = obj.sys.U;
            
            % --- define optimization variables ---
            x = sdpvar(n, N+1);
            z = sdpvar(n, N+1);
            x_0 = sdpvar(n, 1);
            u = sdpvar(m, N);
            
            % define parameter decision variables
            Delta_d=sdpvar(p,1);
            
            % --- define objective ---
            objective = 0;
            for i = 1:N
                objective = objective + x(:,i)'*Q*x(:,i) + u(:,i)'*R*u(:,i);
            end
            objective = objective + x(:,N+1)'*P*x(:,N+1);
            
            % --- define constraints ---
            constraints=[z(:,1)==x_0, x(:,1)==x_0];
            for i = 1:N
                %%% TODO %%%
                % Implement the dynamics which use the oracle using x
                % --- start inserting here ---      
                constraints = [constraints, x(:,i+1) == 0];
                % --- stop inserting here ---
                %%%
                
                constraints = [constraints, z(:,i+1) == A*z(:,i) + B*u(:,i)];
                if i == 1
                    constraints = [constraints, X.A*z(:,i)<=X.b];
                    constraints = [constraints, U.A*u(:,i)<=U.b];
                else
                    Z_i = X - F{i}; Z_i = Z_i.minHRep();
                    constraints = [constraints, Z_i.A*z(:,i)<=Z_i.b];
                    U_i = U - K*F{i}; U_i = U_i.minHRep();
                    constraints = [constraints, U_i.A*u(:,i)<=U_i.b];
                end
            end
            Z_f = X_f - F{N+1}; Z_f = Z_f.minHRep();
            constraints = [constraints, Z_f.A*z(:,N+1)<=Z_f.b];
            
            % --- setup Yalmip solver object ---
            obj.prob = optimizer(constraints, objective, sdpsettings('solver',solver), {x_0, Delta_d}, {u(:,1)});
        end
        
        function F = compute_robust_tightening(obj, K)
            %COMPUTE TIGHTENING Computes the disturbance reachable sets (DRS)
            %   and the corresponding tightenings.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    K = -dlqr(obj.sys.A, obj.sys.B, obj.params.Q, obj.params.R);
                    
                case 2
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%
            
            % look up horizon length
            N = obj.params.N;
            
            % look up system matrices
            A = obj.sys.A;
            B = obj.sys.B;
            
            % disturbance set
            W = obj.sys.W;
            
            % compute disturbance reachable sets (DRS)
            F{1} = Polyhedron(); % initialize empty set
            for i = 1:N
                F_temp = F{i} + ((A+B*K)^(i-1))*W;
                F{i+1} = F_temp.minHRep(); % we need this to reduce the computational complexity
            end
        end
    end
end