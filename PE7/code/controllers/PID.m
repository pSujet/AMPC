%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021, ETH Zurich, {adidier, jsieber}@ethz.ch
%
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2021 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PID < Controller
    %SMMPC SM-MPC Controller Class
    %   Construct and solve RAMPC Problem
    
    properties
        sys
        pos_error_int
        vel_error_int
        prev_pos_error
        prev_vel_error
        Ts
    end
    
    methods
        function obj = PID(sys, params)
            %SMMPC Construct an instance of this class
            %   Construct SM-MPC Class and initialize solver
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            obj@Controller(sys, params);
            
            obj.sys = sys;
            
            % set integral and previous error to zero
            obj.pos_error_int = 0;
            obj.vel_error_int = 0;
            obj.prev_pos_error = 0;
            obj.prev_vel_error = 0;
            
            % get sampling time
            obj.Ts = sys.params.Ts;
        end
        
        function [u, info] = solve(obj, x)
            %COMPUTE_TUBE Computes tube & terminal controller and terminal
            %cost. The method is from Appendix A in Köhler et al. (2019),
            %"Linear robust adaptive model predictive control:
            %Computational complexity and conservatism"
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise
                    error('Wrong number of inputs!')
            end
            
            %%%
            
            %%% Position Error
            % proportional term
            pos_error = obj.params.setpoint(1) - x(1);
            P = pos_error;
            
            % integral term
            obj.pos_error_int = obj.pos_error_int + pos_error*obj.Ts;
            I = obj.pos_error_int;
            
            % derivative term
            D = (pos_error - obj.prev_pos_error)/obj.Ts;
            obj.prev_pos_error = pos_error;
            
            u_pos = obj.params.Kp*P + obj.params.Ki*I + obj.params.Kd*D;
            
            %%% Velocity Error
            % proportional term
            vel_error = obj.params.setpoint(2) - x(2);
            P = vel_error;
            
            % integral term
            obj.vel_error_int = obj.vel_error_int + vel_error*obj.Ts;
            I = obj.vel_error_int;
            
            % derivative term
            D = (vel_error - obj.prev_vel_error)/obj.Ts;
            obj.prev_vel_error = vel_error;
            
            u_vel = obj.params.Kp*P + obj.params.Ki*I + obj.params.Kd*D;
            
            %%% Compute Input
            u = u_pos + u_vel;
            info = '';
            
            % this only works for 1D inputs!
            if u > max(obj.sys.U.V)
                u = max(obj.sys.U.V);
                info = 'Input capped to upper constraint.';
            elseif u < min(obj.sys.U.V)
                u = min(obj.sys.U.V);
                info = 'Input capped to lower constraint.';
            end
        end
    end
end

