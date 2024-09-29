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

classdef BO < handle
    %BAYESIAN OPTIMIZATION BO class implementing a generic Bayesian
    %optimization algorithm.
    properties
        params %parameter struct
        gp %surrogate for objective we model (here GP)
    end
    methods
        function obj = BO(data, params)
            %BO Bayesian Optimization Class Constructor
            %   constructs BO object according to the provided
            %   parameters
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            obj.params = params;
            
            % initialize GP with provided data and kernel
            obj.gp = GP(params, data, params.kernel);
        end
        
        function [theta, y_hat] = acquisition_fn(obj,X)
            %ACQUISITION_FUNCTION Sample acquisition function at provided
            %sampling points X.
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            %%% TODO: Implement LCB as acquisition function
            %%% return the sampled acquisition function in y_hat (this is
            %%% only needed for visualization) and the next parameter to
            %%% sample in theta.
            % --- start inserting code here --- 
            [y, var] = obj.gp.predict(X);
            beta = 1;
            y_hat = y-beta*sqrt(diag(var));
            [~,idx] = min(y_hat);            
            theta = X(:,idx);
            % --- stop inserting code here ---
        end
        
        function add_data(obj,data)
            %ADD_DATA Add new data points to GP
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            obj.gp.add_data(data);
        end
        
        function [theta, mean, conf] = get_estimate(obj,X)
            %GET_ESTIMATE Sample GP at provided sampling points X and
            %return best parameter estimate, expected value of learned
            %function and 95% confidence interval
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            [y, var] = obj.gp.predict(X);
            [mean, idx] = min(y);
            theta = X(:,idx);
            std = sqrt(diag(var));
            conf = [mean-1.96*std(idx), mean+1.96*std(idx)];
        end
        
        function [y, std] = sample(obj,X)
            %SAMPLE Sample GP at sampling points X
            
            %%% Parse inputs %%%
            switch nargin
                case 2
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            [y, var] = obj.gp.predict(X);
            std = sqrt(diag(var));
        end
        
        function data = get_data(obj)
            %GET_DATA Return data stored in GP
            
            %%% Parse inputs %%%
            switch nargin
                case 1
                    
                otherwise 
                    error('Wrong number of inputs!')
            end
            %%%%%%%%%%%%%%%%%%%
            
            data = obj.gp.data;
        end
    end
end