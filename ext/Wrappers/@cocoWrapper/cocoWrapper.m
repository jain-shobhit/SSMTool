classdef cocoWrapper < matlab.mixin.SetGet
    %COCOWRAPPER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        system     % dynamical system object 
        nCycles    % number of periods in (forward) steady state simulation
        outdof     % out dof for result visualization
        periodsRatio = 1  % response period / excitation period
        Options = cocoOptions();
        initialGuess = 'forward' % 'linear'
        branchSwitch = false % true
        multiFnl = []
    end
    
    methods
        % constructor
        function obj = cocoWrapper(sys, ncyles, outdof)
            %COCOWRAPPER Construct an instance of this calss
            obj.system  = sys;
            obj.nCycles = ncyles;
            obj.outdof  = outdof;
        end
        
        function set.initialGuess(obj,initialGuess)
            switch lower(initialGuess)
                case 'forward'
                    obj.initialGuess = 'forward';
                case 'linear'
                    obj.initialGuess = 'linear';
                otherwise
                    error('Unknown omegaSampStyle type: set uniform or nonuniform frequency sampling types')
            end
        end 
        
        % convert fnl from tensor format to multiindex
        function fnlTensor2Multi(obj)
            fnl = obj.system.fnl;
            y = tensor_to_multi_index(fnl);
            obj.multiFnl = y;
        end

        
        % vector field compatiable to coco - autonomous case
        f = aut_ode(obj, x, p, data)
        
        % vector field compatible to coco - nonautonomous case
        f = ode_het(obj, t, x, p, data);
        
        % internal force in M\ddot{u}+N(u,\dot{u}) = F(t,p)
        y = Nhan(obj,u,v);
        
        % derivative of N w.r.t u
        y = dNdu(obj,u,v);
        
        % derivative of N w.r.t v
        y = dNdv(obj,u,v);
        
        % external force in M\ddot{u}+N(u,\dot{u}) = F(t,p)
        y = Fext(obj,t,p);
        
        % derivative of Fext w.r.t p
        y = dFextdp(obj,t,p);
        
        % setup coco options
        prob = cocoSet(obj, prob)
        
        % extract backbone curve for mode whose natural frequency is
        % closest to omega among all natural frequencies
        bd = extract_backbone(obj, omega, varargin) 
        
        % extract Force Response Curve for given frequency range
        bd = extract_FRC(obj, omega_range, varargin)
        
        % extract FRC using forward simulation
        bd = forward_FRC(obj, omega_range, varagin)
    end
end

