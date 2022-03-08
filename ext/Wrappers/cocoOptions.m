classdef cocoOptions < matlab.mixin.SetGet
    %COCOOPTIONS Options for continuatin package coco
    
    properties
        dir_name          % output directory
        
        % settings for continuation
        NPR = 10          % frequency of screen outputs
        NSV = 10          % frequency of storing solutions to disk
        NAdapt = 0        % adaptation period, 0 ==  off
        h0 = 0.1          % initial step size
        h_max = 0.5       % max step size
        h_min = 0.01      % min step size
        h_fac_max = 2     % max step size adaptation factor
        h_fac_min = 0.5   % min step size adaptation factor
        MaxRes = 0.1      % max residual norm in prediction
        bi_direct = true  % go in both directions or not
        PtMX = 100        % max continuation step
        al_max = 7        % max angle between consecutive tangents
        
        % settings for correction
        ItMX = 10         % max. number of iterations
        TOL  = 1e-6       % tolerance for Newton iteration
        
        % settins for collocation
        NTST = 10         % number of mesh intervals
        NCOL = 4          % number of collocation points
        MXCL = true       % enable/disable termination when discretization error exceeds tolerance
    end
    
    methods
        % write your functions here
    end
end

