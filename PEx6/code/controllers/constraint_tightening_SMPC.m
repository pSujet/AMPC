%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2022 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef constraint_tightening_SMPC < Controller
    %ct-SMPC Construct an instance of this class
    %   Construct the stochastic constraint tightening MPC Class
    
    properties
        solver %solver name
        
        % since this controller uses feedback, we also store K
        K % feedback gain
        
        % since we need to compute the RPI set, we additionally store the
        % system information
        sys % System object
    end
    
    methods
        function obj = constraint_tightening_SMPC(sys, params, K, P, p, solver)
            %ct-SMPC Construct an instance of this class
            %   Construct the stochastic constraint tightening MPC Class and
            %   initialize solver
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    [K, P] = dlqr(sys.A, sys.B, params.Q, params.R);
                    K = -K;
                    p = 0.9;
                    solver = 'quadprog';
                    
                case 3
                    error('Wrong number of inputs!')
                    
                case 4
                    p = 0.9;
                    solver = 'quadprog';
                    
                case 5
                    solver = 'quadprog';
                    
                case 6
                    
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
            
            %%% Initialize the constraint tightening SMPC Controller
            % --- compute tightening ---
            obj.K = K;
            % compute robust tightening
            F = obj.compute_robust_tightening(K);
            % compute stochastic tightening
            [Fw_x, Fw_u] = obj.compute_stochastic_tightening(p);
            
            % --- compute terminal set ---
            X_f = obj.compute_MRPI(K);
            
            % look up horizon length
            N = obj.params.N;
            
            % look up state and input dimensions
            n = obj.sys.n;
            m = obj.sys.m;
            
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
            x_bar = sdpvar(n, N+1);
            x_0 = sdpvar(n, 1);
            u_bar = sdpvar(m, N);
            
            % --- define objective ---
            objective = 0;
            for i = 1:N
                objective = objective + x_bar(:,i)'*Q*x_bar(:,i) + u_bar(:,i)'*R*u_bar(:,i);
            end
            objective = objective + x_bar(:,N+1)'*P*x_bar(:,N+1);
                       
            % --- define constraints ---
            constraints = [x_bar(:,1) == x_0];
            for i = 1:N
                constraints = [constraints, x_bar(:,i+1) == A*x_bar(:,i) + B*u_bar(:,i)];
                if i == 1 % in first time step we have no constraints
                    % do nothing
                elseif i == 2 % in second time step we have only stochastic tightening
                    Z_i = X - Fw_x; Z_i = Z_i.minHRep();
                    constraints = [constraints, Z_i.A*x_bar(:,i)<=Z_i.b];
                    V_i = U - Fw_u; V_i = V_i.minHRep();
                    constraints = [constraints, V_i.A*u_bar(:,i)<=V_i.b];
                else
                    Z_i = X - (A+B*K)*F{i-1} - Fw_x; Z_i = Z_i.minHRep();
                    constraints = [constraints, Z_i.A*x_bar(:,i)<=Z_i.b];
                    V_i = U - K*(A+B*K)*F{i-1} - Fw_u; V_i = V_i.minHRep();
                    constraints = [constraints, V_i.A*u_bar(:,i)<=V_i.b];    
                end
            end
            Z_f = X_f - F{N+1}; Z_f = Z_f.minHRep();
            constraints = [constraints, Z_f.A*x_bar(:,N+1)<=Z_f.b];
            
            % --- setup Yalmip solver object ---
            obj.prob=optimizer(constraints, objective, sdpsettings('solver',solver), {x_0}, {u_bar(:,1)});    
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
            
            % look up disturbance set
            W = obj.sys.W;
            
           % compute disturbance reachable sets (DRS)
            F{1} = Polyhedron(); % initialize empty set
            for i = 1:N
                F_temp = F{i} + ((A+B*K)^(i-1))*W;
                F{i+1} = F_temp.minHRep(); % we need this to reduce the computational complexity
            end
        end
        
        function [Fw_x, Fw_u] = compute_stochastic_tightening(obj, p)
            %COMPUTE TIGHTENING Computes the stochastic backoff term for a
            %   uniformly distributed disturbance.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    p = 0.9;
                    
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%
            
            % look up disturbance set
            W = obj.sys.W;
            
            % compute state backoff term
            Fw_x = sqrt(p)*W;
            
            % compute input backoff term
            %temp = obj.K*W;
            Fw_u = obj.K*Fw_x;
        end
        
        function F = compute_MRPI(obj, K)
            %COMPUTE MRPI Computes maximal robust positive invariant set of
            %   the system.
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            MRPI = ULTISystem('A', obj.sys.A + obj.sys.B*K, 'E', eye(obj.sys.n));
            MRPI.x.min = [-obj.sys.X.b(2); -obj.sys.X.b(4)];
            MRPI.x.max = [obj.sys.X.b(1); obj.sys.X.b(3)];
            MRPI.d.min = [-obj.sys.W.b(2); -obj.sys.W.b(4)];
            MRPI.d.max = [obj.sys.W.b(1); obj.sys.W.b(3)];

            disp("Computing MRPI set...");
            F = MRPI.invariantSet();
            disp("done!");
        end
    end
end

