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

function plot_PE8(figNr, nrSteps, X, U, P, gamma, params)
%PLOT PE 8 Summary of this function goes here
%   Detailed explanation goes here

    %%% Parse input arguments %%%
    switch nargin
        case 2
            X = [];
            U = [];
            P = [];
            gamma = [];
            params = {};

        case 3
            U = [];
            P = [];
            gamma = [];
            params = {};
        
        case 4
            P = [];
            gamma = [];
            params = {};
            
        case 5
            gamma = 1;
            params = {};
            
        case 6
            params = {};
            
        case 7
            
        otherwise
            error('Wrong number of inputs!')
    end
    %%%%%%%%%%%%%%%%%%%
    
    figure(figNr)
    hold on
    if ~isempty(X)
        X = X*(180/pi);
        p2 = X.plot('wire', true, 'linestyle', '-', 'linewidth', 2);
    else
        p2 = [];
    end
    if ~isempty(P) && ~isempty(gamma)
        if ~iscell(P)
            x=sdpvar(2,1);
            p3 = YSet(x,[x'*P*(pi/180)^2*x<=gamma]).plot('Color','b','alpha',0.3,'grid',100);
        end
    else
        p3 = [];
    end
    xlabel('position [deg]')
    ylabel('velocity [deg/s]')
    legend([p2,p3],{'constraints', 'safe set'}, 'Location','ne')
    grid()
    set(gcf,'position',[100,100,params.width,params.height],'color','white')
end

