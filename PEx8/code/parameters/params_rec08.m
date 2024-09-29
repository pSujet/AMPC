%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2022, ETH Zurich, {adidier, jsieber}@ethz.ch
% 
% This code is only made available for students taking the advanced MPC class
% in the fall semester of 2022 (151-0371-00L) and is NOT to be distributed.
%
% Authors: Alexandre Didier, Jérôme Sieber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = params_rec08()
%% System Paramters
params.sys.generateNoise = @gaussian;
params.sys.noiseMean = 0;
params.sys.noiseCovariance = 2;
params.sys.noiseArgs = {params.sys.noiseMean, params.sys.noiseCovariance};

%% Plot Parameters
params.plot.height = 350;
params.plot.width = 900;
params.plot.alpha = 0.2;
params.plot.color = [0.7216, 0.1490, 0.0039];
params.plot.lw = 1;
end