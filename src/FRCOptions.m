classdef FRCOptions < matlab.mixin.SetGet
    %ManifoldOptions This class sets the different options available for
    %the Manifold class
    
    properties        
        nRho = 100                 % number of discrete rho values in the range [0, rhomax] (relevant for method == 'level set') 
        nPar = 100                 % number of discrete parameter (Omega/epsilon) values in the parameter range 
        nPsi = 100                 % number of discrete psi values in the range [0, 2*pi] (relevant for method == 'level set') 
        nt = 128                   % number of discrete time intervals over which the periodic orbit is discretized (relevant for post-processing only)
        rhoScale = 1               % factor for increasing rhomax polar FRC
        nCycle = 200               % number of time integration cycles (relevant for numerical time integration)
        
        omegaSampStyle = 'uniform' % 'uniform', 'cocoOut', 'cocoBD' (relevant for method == 'continuation ep/po')
        initialSolver = 'forward'  % 'forward', 'fsolve' (relevant for method == 'continuation ep/po')
        coordinates = 'polar'      % coordinates for solving reduced dynamics: 'polar', 'cartesian'
        sampStyle = 'cocoBD'       % 'uniform', 'cocoOut', 'cocoBD' (relevant for method == 'continuation ep/po')
        method = 'level set'       % 'level set', 'continuation ep', 'continuation po' 
        saveIC = true              % whether save initial conditions on periodic orbit or not 
        init                       % initial solution guess (relevant for method == 'continuation ep/po')
        frac   = [1 1]             % [frac1,frac2] can be used to
                                   % tune the range of subintervals. Specifically, [oma, omb] will
                                   % be changed as [frac1*oma, frac2*omb] except on the end points of frequency range.
        outdof = []                % output degree-of-freedom  
        
        p0 = [];                   % parameters (epsilon,omega) in initial solution guess used in continuation
        z0 = [];                   % states (in slow-time reduced dynamics) in initial solution guess used in continuation
        
        nonAutoParRedCom = false   % compute_perturbed_wisker is called in the parallel computation of non-autonomous SSMs for
                                   % each sampled excitation frequency. The manifold object is transferred in such a call. This 
                                   % communication cost is intenstive in parallel computation. To reduce the communication load,
                                   % we calculate the non-autonomous directly, instead of calling the routine if the field is true 
                                   
        parSamps = []              % solution at specific parameter values used for comparison or verification
        torRotDiret = 'pos'        % 'pos', 'neg' (rotation direction of tori)
        torNumSegs  = 10           % number of Fourier modes in the approximation of tori
        torPurtb    = 1e-4         % perturbation to Neimark-Sacker periodic orbits to yield initial tori
    end
    
    methods      
        
        function set.omegaSampStyle(obj,omegaSampStyle)
            switch lower(omegaSampStyle)
                case 'uniform'
                    obj.omegaSampStyle = 'uniform';
                case 'cocoout'
                    obj.omegaSampStyle = 'cocoOut';
                case 'cocobd'
                    obj.omegaSampStyle = 'cocoBD';
                otherwise
                    error('Unknown omegaSampStyle type: set uniform or nonuniform frequency sampling types')
            end
        end
        
        function set.initialSolver(obj,initialSolver)
            switch lower(initialSolver)
                case 'forward'
                    obj.initialSolver = 'forward';
                case 'fsolve'
                    obj.initialSolver = 'fsolve';
                otherwise
                    error('Unknown omegaSampStyle type: set uniform or nonuniform frequency sampling types')
            end
        end 
        
        function set.coordinates(obj,coordinates)
            switch lower(coordinates)
                case 'polar'
                    obj.coordinates = 'polar';
                case 'cartesian'
                    obj.coordinates = 'cartesian';
                otherwise
                    error('Unknown torus calculation method: set forward/interpolation/pde')
            end
        end 
        
        function set.sampStyle(obj,sampStyle)
            switch lower(sampStyle)
                case 'uniform'
                    obj.sampStyle = 'uniform';
                case 'cocoout'
                    obj.sampStyle = 'cocoOut';
                case 'cocobd'
                    obj.sampStyle = 'cocoBD';
                otherwise
                    error('Unknown sampStyle type: set uniform/cocoOut/cocoBD frequency sampling types')
            end
        end 

        
        function set.method(obj,method)
            switch lower(method)
                case 'level set'
                    obj.method = 'level set';
                case 'continuation ep'
                    obj.method = 'continuation ep';
                case 'continuation po'
                    obj.method = 'continuation po';
                otherwise
                    error('Unknown "method": set "level set" / "continuation ep" / "continuation po" types')
            end
        end        
    end
end

