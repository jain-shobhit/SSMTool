classdef SSM < Manifold
    %SSM Spectral Submanifold Computation
    %   Detailed explanation goes here
    
    properties
       contOptions = cocoOptions();
       FRCOptions = FRCOptions();
       FRSOptions = FRSOptions();
       
       R_1 = []; % coefficients nonaut reduced dynamics, used in 'continuation po'
       W_1 = []; % nonautonomoues ssm coefficients, used in 'continuation po'
    end
    
    methods
        %% other methods
           
        FRC = extract_FRC(obj, parName, parRange, order)               
        
        BB = extract_backbone(obj,modes,order,outdof,varargin)
        
        varargout = FRC_cont_ep(obj,oid,modes,order,mFreq,parName,parRange,outdof,varargin);
        
        [FRC] = FRC_level_set(obj, resMode, order, parName, parRange) 

        varargout = FRC_cont_po(obj,oid,resModes,order,parRange);        
        
        [SD] = extract_Stability_Diagram(obj,resModes, order, omRange, epsRange, parName, p0, varargin)
        
        FRCs = SSM_poSweeps(obj,oid,resonant_modes,order,mFreqs,epSamp,omRange,varargin)
        
        FRCs = SSM_lvlSweeps(obj, omRange, epsSamp, ORDER)
        
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
        % SSM-tor toolbox
        SSM_TR2tor(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_tor2tor(obj,oid,run,lab,parName,parRange,outdof,varargin);
        SSM_BP2tor(obj,oid,run,lab,parName,parRange,outdof,varargin);        
        varargout = auto_po_solver(obj,R_0,oid,t0,p0,coordinates);
        % parallel computation
        activate_parallel(obj,varargin);
        deactivate_parallel(obj);
        % extract frequency response surface
        extract_FRS(obj,oid,modes,order,mFreq,parRange,outdof,optdof,scale_state,scale_obs,varargin);
        % extract ridges and trenches along a frequency response surface
        % without computing it
        varargout = extract_ridges_trenches(obj,oid,resonant_modes,order,mFreqs,parRange,outdof,optdof,varargin);        
    end
    
    methods (Access = protected)
        FRC = SSM_cont_ep(obj,type,oid,run,lab,parName,parRange,outdof,varargin);
        varargout = SSM_cont_po(obj,type,oid,run,lab,initsol,parName,parRange,outdof,varargin);
        varargout = SSM_cont_tor(obj,type,oid,run,lab,parName,parRange,outdof,varargin)
    end
end

