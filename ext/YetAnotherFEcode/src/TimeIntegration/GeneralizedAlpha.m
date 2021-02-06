classdef GeneralizedAlpha < handle
%% Class for Generalized-alpha scheme for time Integration

% We use the following reference in this implementation 
% [1] Arnold, M., & Brüls, O. (2007). Convergence of the generalized-? scheme 
% for constrained mechanical systems. Multibody System Dynamics, 18(2), 185–202. 
% https://doi.org/10.1007/s11044-007-9084-0

% For the time being, we consider systems without constraints.
    properties
        rho_inf = 0.7 % 0<spectral radius<1: rho_inf = 1 means undamped, rho_inf = 0 means asymptotic annhilation of high frequency response
        alpha_m
        alpha_f
        beta
        beta_p
        h = 0 % time step size
        gamma
        gamma_p
        tol = 1e-6
        Solution
        MaxNRit = 10
        ATS = false     % where adaptive time stepping should be on (true) or not (false)
        hmin = 0            % minimum timestep size
        NROpt = 3       % Maximum no. of N-R Iterations
        linear = false
    end
    methods
        function TI = GeneralizedAlpha(varargin)
            %% Input parsing
            p = inputParser;
            addParameter(p,'rho_inf',TI.rho_inf, @(x) isscalar(x) && x<1 && x>0 );
            addParameter(p,'timestep',TI.h, @(x)validateattributes(x, ...
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
            % main time integration parameter that determines numerical
            % damping of high frequency modes.
            TI.rho_inf = p.Results.rho_inf;  
            
            % from (24) in Arnold & Bruels [1], optimal parameters deduced
            % by Chung & Hulbert.
            TI.alpha_m = (2 * TI.rho_inf - 1) / (TI.rho_inf + 1);            
            TI.alpha_f =  TI.rho_inf / (TI.rho_inf + 1);            
            % from (20) in Arnold & Bruels [1] to ensure second-order truncation error 
            TI.gamma = 0.5 + TI.alpha_f - TI.alpha_m;
            TI.beta = 0.25 * (TI.gamma + 0.5)^2;
                        
            % Other parameters by user
            TI.h = p.Results.timestep;
            TI.tol = p.Results.RelTol;
            TI.MaxNRit = p.Results.MaxNRit;
            TI.hmin = p.Results.hmin;
            TI.linear = p.Results.linear;
            
            % Auxillary parameters (these definition are different from the
            % ones given in Arnold & Bruels [1]), we solve for acceleration
            % corrections here instead of displacements.
            TI.beta_p = TI.h^2 * TI.beta * (1 - TI.alpha_f)/(1-TI.alpha_m);
            TI.gamma_p = TI.gamma * TI.h * (1 - TI.alpha_f) / (1 - TI.alpha_m);
        end
                                
        function Integrate(obj,x0,xd0,xdd0,tmax,Residual)
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
            a = xdd0;
            
            while t < tmax
                t = t+obj.h;
                i = i+1;
                [q_new,qd_new,qdd_new,a] = obj.Prediction(q_old,qd_old,qdd_old,a);                
                
                it = -1; % iteration counter
                %% linear case
                if obj.linear 
                    it = it + 1; 
                    % from Algorithm 1 in Arnold & Bruels [1]
                    [r, M, C_t, K_t] = Residual(q_new,qd_new,qdd_new,t);
                    S =  M  + obj.gamma_p * C_t + obj. beta_p * K_t;
                    Da = -S\r;
                    [q_new,qd_new,qdd_new] = obj.Correction(q_new,qd_new,qdd_new,Da);
                    epsilon = 0;
                %% Nonlinear case    
                else 
                    %% Newton-Raphson iterations
                    while true                         
                        it = it + 1;
                        
                        %% Compute Residual and Tangent operators
                        [r, M, C_t, K_t, c0] = Residual(q_new,qd_new,qdd_new,t);                        
                        
                        %% Check convergence
                        epsilon = norm(r)/c0;
                        disp(['Iteration ' num2str(it) ', Residual norm = '  num2str(epsilon)])
                        if (epsilon<obj.tol)  % Error < Tolerance : break
                            break;
                        else % Error >= Tolerance : perform correction from Algorithm 1 in Arnold & Bruels [1]
                            S =  M  + obj.gamma_p * C_t + obj.beta_p * K_t;
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
                a = a + qdd_new * (1-obj.alpha_f)/(1 - obj.alpha_m);

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
        
        function[q,qd,qdd,a] = Prediction(obj,q0,qd0,qdd0,a0)
            % from Algorithm 1 in Arnold & Bruels [1]
            qd = qd0 + obj.h * (1 - obj.gamma) * a0;
            q = q0 + obj.h * qd0 + (0.5-obj.beta) * obj.h^2 * a0;
            a = 1/(1-obj.alpha_m)*(obj.alpha_f * qdd0 - obj.alpha_m * a0);
            q = q + obj.h^2 * obj.beta * a;
            qd = qd + obj.h * obj.gamma * a;  
            qdd = zeros(length(q0),1);
        end
        
        function [q,qd,qdd] = Correction(obj,q,qd,qdd,Da)
            % from Algorithm 1 in Arnold & Bruels [1]
            q = q + obj.beta_p * Da;
            qd = qd + obj.gamma_p * Da;
            qdd = qdd + Da;
        end
    end
end
