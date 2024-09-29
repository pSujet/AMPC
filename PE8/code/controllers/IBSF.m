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

classdef IBSF 
    %IBSF Safety Filter Class
    %   Class for safety filter 
    
    properties
        params
        sys
        K
        P
    end
    
    methods
        function obj = IBSF(sys, params, solver)
            %IBSF Construct an instance of this class
            %   Construct safety filter Class and initialize solver
            %%% Parse inputs %%%
            switch nargin
                case 2
                    solver = '';
                    
                case 3
                    
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            obj.params = params;
            
            obj.sys = sys;
            
            %%% Execise 2: Invariant Set Computations
            E = sdpvar(sys.n, sys.n);
            Y = sdpvar(sys.m, sys.n);
            
            % Lookup table
            A = sys.A;
            B = sys.B;
            X = sys.X;
            U = sys.U;
            
            %%% TODO: Implement an SDP to compute the largest ellipsoidal
            %%% invariant set which satisfies the state constraints and the 
            %%% corresponding safe controller which satisfies the input 
            %%% constraints. You can use logdet(.) in the objective and the 
            %%% invariance, state and input constraint formulations seen in 
            %%% the lecture.
            % --- start inserting code here ---
            
            objective = -log(det(E));
            constraints = [E >= 1e-3*eye(sys.n)];
            constraints = [constraints, [E, (E*A'+Y'*B'); (A*E + B*Y), E]>=0];
            % state constraints
            for i = 1:size(sys.X.b,1)
                constraints = [constraints,sys.X.A(i,:)*E*sys.X.A(i,:)' <= sys.X.b(i)^2];
            end
            % input constraints
            for i = 1:size(sys.U.b,1)
                constraints = [constraints,[sys.U.b(i)^2, sys.U.A(i,:)*Y;Y'*sys.U.A(i,:)' E] >= 0];
            end                
           
            optimize(constraints, objective, sdpsettings('verbose',0,'solver',solver));
            
            obj.P = inv(value(E));
            obj.K = value(Y)*obj.P;
            % --- stop inserting code here ---
            
        end
        
        function [u_out] = solve(obj, x, u)
            %SOLVE Returns input to be applied
            % Performs the safety check and returns the input to be applied
            
             %%% Parse inputs %%%
            switch nargin
                case 3
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            % Exercise 3 and 4 
            
            %%% TODO: Implement the safety filter checks provided in the
            %%% exercise sheet. The inputs to this function are x, which is
            %%% the current state x(k) and u, which is the proposed
            %%% learning-based input.
            % --- start inserting code here ---
            A = obj.sys.A;
            B = obj.sys.B;
            if (A*x+B*u)'*obj.P*(A*x+B*u)<= 1
                u_out = u;
            else
                u_out = obj.K*x;
            end
            
            
            % --- stop inserting code here ---
        end
    end
end

