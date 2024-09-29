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

classdef SM < Parameter_Estimator
    %SM Set-membership Estimation Class
    % Class used to compute the set-membership update
    
    properties
        Omega % set of parameters
        theta_hat % center of parameter set estimate
        sys
    end
    
    methods
        function obj = SM(sys)
            %SM Construct an instance of this class
            %   Set-membership estimation class
            %   Inputs: System object and bool if sets are hypercubes
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            obj@Parameter_Estimator(sys);
            
            obj.sys = sys;
            
            obj.Omega = Polyhedron(obj.params.A_theta, obj.params.b_theta);
            obj.theta_hat = sys.theta;
        end
        
        function [Delta]=update_estimate(obj, x, x_prev, u_prev)
            %UPDATE_ESTIMATE Compute Parameter Set Estimate at the next
            %time step
            %   Compute Parameter Set Estimate at the next
            %time step given a new measurement x, a previous measurement
            %x_prev, the applied input and the system
            
            %%% Parse inputs %%%
            switch nargin
                case 4
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            %%% TODO
            % Implement the set-membership update. Please note that you also
            % need to implement the redundant constraint removal method
            % before you run the set-membership estimator for the first
            % time.
            
            % --------- Start Modifying Code Here -----------            
            A_Omega = obj.params.A_theta;
            b_Omega = obj.params.b_theta;
            
            % compute non-falsified parameter set
            A_Delta = - obj.params.A_w * (obj.params.A_delta{2} * x_prev(1) + obj.params.A_delta{3}*x_prev(2));
            b_Delta = obj.params.b_w - obj.params.A_w*(x - obj.params.A_delta{1}*x_prev - obj.params.B_delta{1}*u_prev);
            
            Delta = Polyhedron(A_Delta, b_Delta);
            
            % compute intersection
            A_Omega = [A_Omega;A_Delta];
            b_Omega = [b_Omega;b_Delta];
            % next also implement this method below!
            [A_Omega, b_Omega] = obj.remove_redundant_halfspace(A_Omega, b_Omega);
            
            % --------- Stop Modifying Code Here -----------
            
            obj.Omega = Polyhedron(A_Omega, b_Omega);
            obj.theta_hat = obj.point_estimate();
        end
        
        function [A_Omega, b_Omega] = remove_redundant_halfspace(obj, A_Omega, b_Omega)
            %REMOVE_REDUNDANT_HALFSPACE Removes redundant halfspaces from
            %A_Omega and b_Omega
            % Returns minimal realization halfspace matrix and vector which
            % describe Omega
            
            %%% Parse inputs %%%
            switch nargin
                case 3
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            %%% TODO
            % Remove redundant halfspaces. You can either use the method
            % presented in the lecture/recitation or suitable MPT3
            % functionalities.
            % Note: If you use the method presented in the lecture:
            % Make sure to exclude halfspace constraints that where already
            % marked as redundant when solving the linear program.
            % Additionally, check if the linear program attains a solution
            % in each iteration (the objective can grow unbounded).
            
            % --------- Start Modifying Code Here -----------
            P = Polyhedron('A', A_Omega, 'b', b_Omega);            
            P.minHRep();                      
            % change these lines accordingly
            A_Omega=P.A;
            b_Omega=P.b;
            
            % --------- Stop Modifying Code Here -----------
        end
        
        function estimate = point_estimate(obj)
            %POINT_ESTIMATE Compute point estimate in parameter set
            %   Compute point estimate of parameter via center of the outer
            %   approximation of the parameter set estimate.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            outBox = obj.Omega.outerApprox;
            estimate = sum(outBox.A.*outBox.b)'/2;
        end
        
        function W_theta = estimate_W_theta(obj)
            %ESTIMATE_W_THETA Compute point estimate in parameter set
            %   Compute point estimate of parameter via center of the outer
            %   approximation of the parameter set estimate.
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            x_vertices = obj.sys.X.V; %vertices
            W_theta = Polyhedron(zeros(0,2),[]);
            for i=1:size(x_vertices,1)
                phi_x = [obj.sys.A_delta{2}*x_vertices(i,:)', obj.sys.A_delta{3}*x_vertices(i,:)'];
                aux = phi_x*(obj.Omega - obj.sys.theta); % use initial parameter estimate
                W_theta = Polyhedron('V',[W_theta.V; aux.V]).minVRep;
            end
        end
    end
end

