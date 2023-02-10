classdef Manifold < matlab.mixin.SetGetExactNames
    %MANIFOLD Invariant manifold computation around
    %   Detailed explanation goes here
    
    
    properties
        System      % dynamical system object
        dimSystem   % Phase space dimensionality

        E = []      % invariant subspace of the linearized system (data structure)
        N = []      % Invariant torus of the linearized system
        
        resonance   % (near) resonances data structure
        dimManifold % manifold dimensionality
                
        Options = ManifoldOptions()                
    end
    
    properties (SetAccess = private)        
        solInfo = struct('memoryEstimate', [], 'timeEstimate', []) 
        % This data structure stores the solution information:
        % memory consumption estimate in MB at each order
        % computational time estimate in seconds at each order
    end

    
    methods
        %% Constructor
        function obj = Manifold(Sys)
            %MANIFOLD Construct an instance of this class
            %   Detailed explanation goes here
            obj.System = Sys;        
        end        

        %% SET methods

        
        %% GET methods        
        function N = get.dimSystem(obj)
            N = obj.System.N;
        end
        
        function M = get.dimManifold(obj)
            M = 0;
            if ~isempty(obj.E)
                M = M + size(obj.E.basis,2);
            end
            
            if ~isempty(obj.N)
                M = M + numel(obj.N.Omega);
            end
        end
        
        %% other methods
        choose_E(obj,modes) % choose master modal space based on modes
        
        [W_0, R_0] = compute_whisker(obj, order)
        
        [W_0j, R_0j, multi_input] = cohomological_solution(obj, i,  W_0, R_0, multi_input)
        
        [W, f] = compute_perturbed_whisker(obj, order)            
           
        [rho] = compute_analyticity_domain(obj,appr_order) %compute analyticity domain at approximation order appr_order.
        
        err = compute_auto_invariance_error(obj,W,R,rhosamp,orders,ntheta,varargin);
        
        res = compuate_invariance_residual(obj,W0,R0,p,type,varargin)
    end
end

