function varargout = time_integration_transient(obj,Omega,varargin)
% TIME_INTEGRATION_PERIODIC This function is used to perform direct 
% numerical time integration of the dynamical system over a 
% range of forcing frequencies such that transients decay to obtain a
% steady-state periodic response 
% obj: object of DynamicalSystem class
% omegaRange: forcing frequency range
% optional arguments:
% 'nCycles': maximum number of forcing cycles over which transient response is
%           computed to achieve steady-state
% 'nOmega': number of equally spaced forcing frequencies in omegaRange over
%           which the periodic response is computed
% 'integrationMethod':  method to be chosen for time integration: 'ode45'
%                       or 'ode15s' or 'Newmark' or 'Galpha'
% 'outdof': the degree of freedom for which the converged periodic orbit 
%           should be plotted
% 'init': initial condition vector to be used for the first time
%           integration

N = obj.N;
[integrationMethod, nCycles, zInit, outdof, nSteps, Plot, samp] =...
    parse_inputs(N, varargin{:});

ndof = numel(outdof);

if abs(Omega)<1e-10
    T = 2*pi*nCycles;
else
    T = 2*pi/Omega*nCycles;
end
obj.Omega = Omega;
odefun = @(t,x) obj.odefun(t,x);    
residual = @(q,qd,qdd,t)obj.residual(q,qd,qdd,t);

% forward simulation
switch integrationMethod
    case 'ode45'
        % transient time integration to converging to steady-state
%         odeoptions = odeset('RelTol',1e-9,'AbsTol',1e-9);
        [t0, x0] = ode45(odefun, linspace(0,T,nSteps*nCycles+1), zInit);

    case 'ode15s'
%         odeoptions = odeset('RelTol',1e-9,'AbsTol',1e-9);
        [t0, x0] = ode15s(odefun, linspace(0,T,nSteps*nCycles+1), zInit);

    case {'Newmark', 'Galpha'}
        % Instantiate object for nonlinear time integration
        if abs(Omega)<1e-10
            h = 2*pi/nSteps;
        else
            h = 2*pi/Omega/nSteps;
        end
        n = obj.n;
        if strcmp(integrationMethod,'Newmark')
            TI = ImplicitNewmark('timestep',h,'alpha',0.005);
        elseif strcmp(integrationMethod,'Galpha')
            TI = GeneralizedAlpha('timestep',h,'rho_inf',0.7);
        end
        % Residual evaluation function handle            
        q0 = zInit(1:n); qd0 = zInit(n+1:N); 
        qdd0 = obj.M\(obj.compute_fext(0) - obj.C*qd0 - obj.K*q0 ...
            - obj.compute_fnl(q0,qd0));

        % Nonlinear Time Integration
        TI.Integrate(q0,qd0,qdd0,T,residual);

        x0 = full([TI.Solution.q; TI.Solution.qd]).';
        t0 = TI.Solution.time;

    otherwise
        error('Unknown integration method chosen \n Valid choices are ode45, ode15s, Newmark and Galpha')            

end

% save output
save('traj.mat','x0', 't0','Omega')

varargout{1} = t0(1:samp:end);
if outdof(1)>0
    varargout{2} = x0(1:samp:end,outdof);
    varargout{3} = x0(end,:); % final state
end

%% Plot initial PO (check)
if outdof(1)>0 && Plot
    for i=1:ndof
        figure;
        i1 = outdof(i);
        plot(t0, x0(:,i1), 'LineWidth', 2);
        xlabel('time','interpreter','latex'); ylabel('z','interpreter','latex')
        title(['Periodic orbit \Omega= ' num2str(Omega)])
        set(gca, 'LineWidth', 1.5);
        set(gca, 'FontSize', 14);
        axis square        
        saveas(gcf,['po_' num2str(i1), '.fig'])
        close(gcf)
    end
end

end


function [integrationMethod,nCycles,zInit,outdof,nSteps,Plot,samp] = parse_inputs(N,varargin)
%% parsing inputs
defaultinit = zeros(N,1);
defaultncycles = 5000;
defaultIntegrationMethod = 'ode45';
defaultoutdof = 0;
defaultnSteps = 30;
defaultPlot = false;
defaultSampGap = 1;

p = inputParser;
addParameter(p,'nSteps',defaultnSteps, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'outdof',defaultoutdof, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','nonnegative'}) );
addParameter(p,'nCycles',defaultncycles, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'integrationMethod',defaultIntegrationMethod);
addParameter(p,'init',defaultinit);
addParameter(p,'plot',defaultPlot);
addParameter(p,'samp',defaultSampGap);


parse(p,varargin{:});

integrationMethod = p.Results.integrationMethod;
nCycles = p.Results.nCycles;
zInit = p.Results.init;
outdof = p.Results.outdof;
nSteps = p.Results.nSteps;
Plot = p.Results.plot;
samp = p.Results.samp;
end
