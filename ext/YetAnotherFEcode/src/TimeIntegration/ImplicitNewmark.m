classdef ImplicitNewmark < handle
    %% Class for Implicit Newmark time Integration 
    % cf. flowchart 7.21 in the following reference:
    % Geradin, M., & Rixen, D. (2015). Mechanical Vibrations : Theory and Application to Structural Dynamics (Third Edition). 
    properties
        % Time integration parameters
        alpha = 0.005
        beta
        h = 0 % time step size
        gamma
        tol = 1e-6      % Relative error tolerance
        
        
        Solution        % Solution data structure
        MaxNRit = 10
        ATS = false     % where adaptive time stepping should be on (true) or not (false)
        hmin = 0        % minimum timestep size (only used when ATS = true)
        NROpt = 3       % Maximum no. of N-R Iterations
        linear = false  % whether system is linear or not
    end
    
    methods
        function TI = ImplicitNewmark(varargin)
            %% Input parsing
            p = inputParser;
            addParameter(p,'timestep',TI.h, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'alpha',TI.alpha, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'RelTol',TI.tol, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','positive'}) );
            addParameter(p,'MaxNRit',TI.MaxNRit, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty','integer','positive'}) );
            addParameter(p,'linear', TI.linear, @(x)validateattributes(x,{'logical'},{'nonempty'}));
            addParameter(p,'hmin', TI.hmin, @(x)validateattributes(x, ...
                {'numeric'},{'nonempty'}) );
            addParameter(p,'ATS', TI.ATS, @(x)validateattributes(x,{'logical'},{'nonempty'}));
            
            parse(p,varargin{:});
            
            %% Properties assignment
            TI.alpha = p.Results.alpha;
            TI.beta = (1 + TI.alpha)^2/4;
            TI.gamma = 1/2 + TI.alpha;
            
            TI.h = p.Results.timestep;
            TI.tol = p.Results.RelTol;
            TI.MaxNRit = p.Results.MaxNRit;
            TI.hmin = p.Results.hmin;
            TI.linear = p.Results.linear;
        end
        function Integrate(obj,x0,xd0,xdd0,tmax, Residual)            
            % Integrates with Initial condition x0,xd0 from [0 tmax]
            % Residual is a function handle that has the following syntax
            % 
            if obj.h ==0
                error('Please specify a positive time step')
            end
            
            tic
            t=0;
            time = t;
            q = x0;
            qd = xd0;
            q_old = x0;
            qd_old = xd0;
            qdd_old = xdd0;
            NR = 0;
            R = 0;
            i = 1;            
            
            while t < tmax
                t = t+obj.h;
                i = i+1;
                [q_new,qd_new,qdd_new] = obj.Prediction(q_old,qd_old,qdd_old);                
                
                it = -1; % iteration counter
                %% linear case
                if obj.linear 
                    it = it + 1; 
                    [r, drdqdd, drdqd, drdq] = Residual(q_new,qd_new,qdd_new,t);
                    S = drdqdd + obj.gamma * obj.h * drdqd + obj.beta * obj.h^2 * drdq;
                    Da = -S\r;
                    [q_new,qd_new,qdd_new] = obj.Correction(q_new,qd_new,qdd_new,Da);
                    epsilon = 0;
                %% Nonlinear case    
                else 
                    %% Newton-Raphson iterations
                    while true                         
                        it = it + 1;
                        
                        %% Compute Residual and Tangent operators
                        [r, drdqdd, drdqd, drdq, c0] = Residual(q_new,qd_new,qdd_new,t);                        
                        
                        %% Check convergence
                        epsilon = norm(r)/c0;
                        disp(['Iteration ' num2str(it) ', Residual norm = '  num2str(epsilon)])
                        if (epsilon<obj.tol)  % Error < Tolerance : break
                            break;
                        else % Error >= Tolerance : perform correction
                            S = drdqdd + obj.gamma * obj.h * drdqd + obj.beta * obj.h^2 * drdq;
                            Da = -S\r;
                            [q_new,qd_new,qdd_new] = obj.Correction(q_new,qd_new,qdd_new,Da);
                        end
                        
                        %% Adapt time step to maintain an optimal number (obj.NROpt) of N-R iterations 
                        if obj.h > obj.hmin && obj.ATS
                            obj.h = max(obj.h*obj.NROpt/it, obj.hmin);
                        end
                        
                        %% When too many iterations
                        if (it > obj.MaxNRit)
                            warning('Max N-R iterations reached')                            
                            if  epsilon > 1 
                                disp('Exiting time integration: Too high a residual')
                                soltime=toc;
                                obj.Solution.time = time;
                                obj.Solution.q = q;
                                obj.Solution.qd = qd;
                                obj.Solution.NR = NR;
                                obj.Solution.R = R;
                                obj.Solution.soltime = soltime;
                                return
                            else
                                disp('Continuing iterations anyway since the residual is low')                            
                            end
                        end                    

                    end
                end
                
                %% Update solution
                time = [time t];
                NR = [NR it];
                R = [R epsilon];
                disp(['time integration completed: ', num2str(100* t/tmax), '%'])
                
                q = [q q_new];
                qd = [qd qd_new];
                q_old = q_new;
                qd_old = qd_new;
                qdd_old = qdd_new;
                
            end
            soltime = toc;
            obj.Solution.time = time;
            obj.Solution.q = q;
            obj.Solution.qd = qd;
            obj.Solution.NR = NR;
            obj.Solution.R = R;
            obj.Solution.soltime = soltime;
        end
        function[q,qd,qdd] = Prediction(obj,q0,qd0,qdd0)
            qd = qd0 + obj.h * (1 - obj.gamma) * qdd0;
            q = q0 + obj.h * qd0 + (0.5-obj.beta) * obj.h^2 * qdd0;
            qdd = zeros(length(q0),1);
        end
        function [q,qd,qdd] = Correction(obj,q,qd,qdd,Da)
            q = q + obj.beta * obj.h^2 * Da;
            qd = qd + obj.gamma * obj.h * Da;
            qdd = qdd + Da;
        end
    end
end
