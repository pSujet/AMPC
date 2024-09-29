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

classdef LinearAffineSystem < System
    %SYSTEM System class implementing a segway system
    %   The System class defines the interfaces, which will be used by the
    %   controller and parameter classes implemented in the recitations.
    properties
        A_delta %dynamic matrices of affine dynamics parameterization; A = A_delta{0} + \sum A_delta{i} * Omega(i)
        B_delta %dynamic matrices of affine dynamics parameterization; B = B_delta{0} + \sum B_delta{i} * Omega(i)
        Omega %uncertain parameter vector
    end
    methods
        function obj = LinearAffineSystem(params)
            %SYSTEM System Class Constructor
            %   constructs system object according to the provided
            %   parameters
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            % call parent class constructor
            obj@System(params);
            
            % define parameter set
            if ~isempty(obj.params.A_theta) && ~isempty(obj.params.b_theta)
                obj.Omega = Polyhedron(obj.params.A_theta, obj.params.b_theta);
            end

            % define dynamic uncertainty in A & B
            obj.A_delta = obj.params.A_delta;
            obj.B_delta = obj.params.B_delta;
            
            % define system matrices
            obj.A = obj.params.A;
            obj.B = obj.params.B;
            obj.C = obj.params.C;
            obj.D = obj.params.D;

            % define update law
            obj.f = @(x,u) obj.A*x + obj.B*u + obj.generateNoise(obj.noiseArgs{:});
            obj.h = @(x,u) obj.C*x + obj.D*u;
        end
        
        function update_params(obj, params)
            %UPDATE_PARAMS
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
           
            % call parent class constructor
            update_params@System(params);
            
            % define parameter set
            if ~isempty(obj.params.A_theta) && ~isempty(obj.params.b_theta)
                obj.Omega = Polyhedron(obj.params.A_theta, obj.params.b_theta);
            end

            % define dynamic uncertainty in A & B
            obj.A_delta = obj.params.A_delta;
            obj.B_delta = obj.params.B_delta;
            
            % define system matrices
            obj.A = obj.params.A;
            obj.B = obj.params.B;
            obj.C = obj.params.C;
            obj.D = obj.params.D;

            % define update law
            obj.f = @(x,u) obj.A*x + obj.B*u + obj.generateNoise(obj.noiseArgs{:});
            obj.h = @(x,u) obj.C*x + obj.D*u;
        end
    end
end