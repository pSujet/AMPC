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

classdef RMPC < Controller
    %RMPC Robust MPC Controller Class
    %   Construct and solve robust linear MPC Problem
    
    properties
        solver % solver name
        
        % since this controller uses feedback, we also store K
        K % feedback gain
        
        sys %System description
    end
    
    methods
        function obj = RMPC(sys, params, rho, solver)
            %   RMPC Construct an instance of this class
            %   Construct the robust MPC Class and
            %   initialize solver
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    rho = 0.9;
                    solver = 'sedumi';
                    
                case 3
                    solver = 'sedumi';
                    
                case 4
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            % call super class constructor
            obj@Controller(sys, params);
            
            % Internal system description
            obj.sys = sys;
            
            % initialize solver name
            obj.solver = solver;
            
            %%% Initialize the RMPC Controller
            
            % --- define optimization variables ---
            z = sdpvar(sys.n, obj.params.N+1); % Nominal prediction states
            x_k = sdpvar(sys.n, 1); % Current true state, which is the input to the obj.prob optimizer object
            v = sdpvar(sys.m, obj.params.N);   % Nominal prediction inputs
            
            %%% TODO %%%
            % Compute the required tightenings given rho
            % You have to define obj.K!
            % Define the robust MPC problem constraints and objective.
            
            % --- start inserting here ---   
            % compute tightenings
            [x_tight, u_tight, P, K, delta] = compute_tightening(obj, rho);
            obj.K = K;
            
            % extract data
            N = obj.params.N;
            A_x = obj.sys.X.A;
            A_u = obj.sys.U.A;
            nx = size(A_x,1);                     % number of states constraints
            nu = size(A_u,1);          
            
            % Define objective            
            objective = 0;
            for i = 1:N
                objective = objective + z(:,i)'*obj.params.Q*z(:,i) + v(:,i)'*obj.params.R*v(:,i);
            end
            % Define constraints
            constraints = [];
            for i = 1:N
                constraints = [constraints, z(:,i+1) == sys.A*z(:,i)+sys.B*v(:,i)];
                constraints = [constraints, sys.X.A*z(:,i)<=sys.X.b - x_tight];
                constraints = [constraints, sys.U.A*v(:,i)<=sys.U.b - u_tight];
            end
            % terminal constraint
            constraints = [constraints, z(:,N+1)==[0;0]];
            constraints = [constraints, (x_k - z(:,1))'*P*(x_k - z(:,1)) <= delta^2];
            % --- stop inserting here ---
            %%%
            
            % --- setup Yalmip solver object ---
            obj.prob = optimizer(constraints, objective, sdpsettings('solver',solver), {x_k}, {v(:,1),z(:,1)});    
        end
        
        function [x_tight, u_tight, P, K, delta] = compute_tightening(obj, rho)
            %COMPUTE TIGHTENING Computes an RPI set and the corresponding
            %   tightening, which minimizes the constraint tightening.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    rho = 0.9;
                    
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%
             
            %--- setup & solve offline optimization problem ---
            
            %%% TODO %%%
            % Define the SDP for computing a contractive ellipsoid
            % Compute the sublevel set delta which renders the ellipsoid
            % RPI as well as the corresponding state and input constraint
            % tightenings. 
            %
            % Make sure to define every output of the function!
            %
            % We recommend weighing the state tightenings in the objective
            % function by a factor of 50, i.e., 50*sum(c_{x,j}^2) in order
            % to obtain good results. However, any RPI set which you
            % compute with the provided constraints in 1a) and which 
            % results in a recursively feasible RMPC scheme, as well as 
            % feasible for the provided initial condition will be accepted.
            % --- start inserting here ---   
            
            
            % Extract data
            A = obj.sys.A;
            B = obj.sys.B;
            A_x = obj.sys.X.A;
            A_u = obj.sys.U.A;
            W = obj.sys.W;                        % noise polytopes
            nx = size(A_x,1);                     % number of states constraints
            nu = size(A_u,1);                     % number of input constraints
            n = obj.sys.n;                        % number of states
            
            % Define decision variable
            E = sdpvar(n,n);                    % square matrix of states
            Y = sdpvar(1,n);                    % Y = KE
            c_xj_squre = sdpvar(nx,1);          % 
            c_uj_squre = sdpvar(nu,1);          % 
            w_bar_square = sdpvar(1,1);           % squre(max(sqrt(w'*inv(E)*w)))
            
            % Define objective
            objectives = (nx + nu)*w_bar_square;
            for i = 1:nx
                objectives = objectives + c_xj_squre(i,1);
            end
            
            for i = 1:nu
                objectives = objectives + c_uj_squre(i,1);
            end
            
            objectives = (1./(2*(1-rho))).*objectives;
            
            % Define constraints            
            constraints = [E >= eye(n)];
            constraints = [constraints, [rho^2*E, (A*E + B*Y)';A*E + B*Y, E]>=0];
            for i = 1:nx
                constraints = [constraints, [c_xj_squre(i,1), A_x(i,:)*E; E'*A_x(i,:)',E] >= 0];
            end
            for i = 1:nu
                constraints = [constraints, [c_uj_squre(i,1), A_u(i,:)*Y; Y'*A_u(i,:)',E] >= 0];
            end
            for i = 1:size(W.V,1)
                % v_w = W.V(1,:)
                constraints = [constraints, [w_bar_square, W.V(1,:); W.V(1,:)', E] >= 0];
            end
            
            
%             options = sdpsettings('solver','sedumi', 'verbose', 0);
            % Only uncomment if you have Mosek installed, it will result in
            % a nice speedup
            options = sdpsettings('solver','mosek', 'verbose', 0);
            
            % run optimization
            optimize(constraints,objectives);
            
            % get data back
            w_bar = sqrt(value(w_bar_square));
            x_tight = sqrt(value(c_xj_squre))*w_bar*(1/(1-rho));
            u_tight = sqrt(value(c_uj_squre))*w_bar*(1/(1-rho));
            P = inv(value(E));
            K = value(Y)*P;
            delta = w_bar*(1/(1-rho));
            
            
            
            % --- stop inserting here ---
            %%%
            
        end
    end
end

