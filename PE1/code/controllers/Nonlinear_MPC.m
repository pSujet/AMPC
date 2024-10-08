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

classdef Nonlinear_MPC < Controller
    %NLMPC Nonlinear MPC Controller Class
    %   Construct and solve nominal nonlinear MPC Problem
    
    properties
        % we need to define the following properties, for Ipopt to work
        % properly as a class object.
        x %state (optimization variable)
        x_0 %initial state (parameter)
        u %input (optimization variable)
        solver %NLP solver name
        objective %nominal nonlinear MPC objective function
    end
    
    methods
        function obj = Nonlinear_MPC(sys, params, solver)
            %NLMPC Construct an instance of this class
            %   Construct non-linear MPC Class and initialize solver
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    solver = 'ipopt';
                    
                case 3
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            % call super class constructor
            obj@Controller(sys, params);
            
            % initialize NLP solver name
            obj.solver = solver;
            
            %%% TODO %%%
            % Define the nonlinear nominal MPC problem below.
            % Use the Opti stack from CasADi to formulate the optimization
            % problem.
            % You have to define obj.prob!
            
            % --- start inserting here ---
            % initialize Opti stack
            obj.prob = casadi.Opti();
            
            % define optimization variables
            obj.x = obj.prob.variable(sys.n,obj.params.N+1);
            obj.x_0 = obj.prob.parameter(sys.n,1);
            obj.u = obj.prob.variable(sys.m,obj.params.N);           
            
            % define objective
            obj.objective = 0;
            for i = 1:obj.params.N
                obj.objective = obj.objective + obj.x(:,i)'*obj.params.Q*obj.x(:,i) + obj.u(:,i)'*obj.params.R*obj.u(:,i);
            end
            
            % minimize objective 
            obj.prob.minimize(obj.objective)
            
            % define constraints
            obj.prob.subject_to(obj.x(:,1) == obj.x_0); % initial condition
            for i = 1:obj.params.N
                obj.prob.subject_to(obj.x(:,i+1) == sys.f(obj.x(:,i),obj.u(:,i)))
                obj.prob.subject_to(sys.X.A*obj.x(:,i) <= sys.X.b)
                obj.prob.subject_to(sys.U.A*obj.u(:,i) <= sys.U.b)            
            end
            obj.prob.subject_to(obj.x(:,obj.params.N+1)== [0;0]) % terminal condition, == 0? 
                        
            % --- stop inserting here ---
            %%%
            
            % --- setup NLP ---
            % initialize non-verbose
            obj.prob.solver(solver, struct('print_time', 0), struct('print_level', 0));          
        end
        
    end
end


