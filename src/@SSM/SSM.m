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
        
        BB = extract_backbone(obj,modes,order,outdof,varargin)
        
        varargout = FRC_cont_ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        
        [FRC] = FRC_level_set(obj, resMode, order, parName, parRange) 
        % SSM-ep toolbox
        varargout = SSM_isol2ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        varargout = SSM_ep2ep(obj,oid,run,lab,parName,parRange,outdof,varargin);
        varargout = SSM_BP2ep(obj,oid,run,lab,parName,parRange,outdof,varargin);
        varargout = SSM_ep2SN(obj,oid,run,lab,parRange,outdof,varargin);
        varargout = SSM_ep2HB(obj,oid,run,lab,parRange,outdof,varargin);
        varargout = SSM_epSweeps(obj,oid,run,lab,epSamps,omRange,outdof,varargin);
        % SSM-po toolbox
        SSM_isol2po(obj,oid,run,lab,initsol,parName,parRange,outdof,varargin);
        SSM_HB2po(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_po2po(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_BP2po(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_po2TR(obj,oid,run,lab,parRange,outdof,varargin);
        SSM_po2SN(obj,oid,run,lab,parRange,outdof,varargin);
        SSM_po2PD(obj,oid,run,lab,parRange,outdof,varargin);
        SSM_po2Tinf(obj,oid,run,lab,parRange,outdof,varargin);
        % SSM-tor toolbox
        SSM_TR2tor(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_tor2tor(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_BP2tor(obj,oid,run,lab,parName,parRange,outdof,varargin);        
        varargout = auto_po_solver(obj,R_0,oid,t0,p0,coordinates);
        % parallel computation
        activate_parallel(obj,varargin);
        deactivate_parallel(obj);
    end
    
    methods (Access = protected)
        FRC = SSM_cont_ep(obj,type,oid,run,lab,parName,parRange,outdof,varargin);
        varargout = SSM_cont_po(obj,type,oid,run,lab,initsol,parName,parRange,outdof,varargin);
        varargout = SSM_cont_tor(obj,type,oid,run,lab,parName,parRange,outdof,varargin)
    end
end

