classdef SSM < Manifold
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
       contOptions = cocoOptions();
       FRCOptions = FRCOptions();
    end
    
    methods
        %% other methods
           
        FRC = extract_FRC(obj, parName, parRange, order)               
        
        BB = extract_backbone(obj,modes,order,outdof)
        
        varargout = FRC_cont_ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        
        [FRC] = FRC_level_set(obj, resMode, order, parName, parRange) 
        
        varargout = SSM_isol2ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        varargout = SSM_ep2ep(obj,oid,run,lab,parName,parRange,outdof,varargin);
        varargout = SSM_BP2ep(obj,oid,run,lab,parName,parRange,outdof,varargin);
        varargout = SSM_ep2SN(obj,oid,run,lab,parRange,outdof,varargin);
        varargout = SSM_ep2HB(obj,oid,run,lab,parRange,outdof,varargin);
        varargout = SSM_epSweeps(obj,oid,run,lab,epSamps,omRange,outdof,varargin);
    end
end

