function [p0,z0] = initial_fixed_point(p0,initialSolver,ispolar,odefun,nCycle,m,varargin)
% INITIAL_FIXED_POINT This function construct initial solution to the fixed
% point of leading-order reduced dynamics. Two methods: forward simulation
% and optimization, are avaliable to obtain such an initial fixed point.
%

if numel(varargin)>0 && iscell(varargin{1})
    p0 = varargin{1}{1};
    z0 = varargin{1}{2};
    z0 = z0(:);
else
    if ispolar
        z0 = 0.1*ones(2*m,1);
    else
        z0 = zeros(2*m,1);
        % solving linear equations
%         for i=1:numel(iNonauto)
%             id  = iNonauto(i);
%             r   = rNonauto(i);
%             rRe = real(r);
%             rIm = imag(r);
%             ai  = fdata.lamdRe(id);
%             bi  = fdata.lamdIm(id)-mFreqs(id)*p0(1);
%             z0i = [ai -bi;bi ai]\[-rRe;-rIm];
%             z0(2*id-1:2*id) = z0i;
%         end
    end
end
% construct initial guess equilibrium points
switch initialSolver
    case 'fsolve'
        % fsolve to approximate equilibrium
        fsolveOptions = optimoptions('fsolve','MaxFunctionEvaluations',100000,...
            'MaxIterations',1000000,'FunctionTolerance', 1e-10,...
            'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-10);
        z0 = fsolve(@(z) odefun(z,p0),z0,fsolveOptions);
    case 'forward'
        % forward simulation to approach equilibirum
        tspan = [0 nCycle*2*pi/p0(1)]; %nCycle
        odefw = @(t,z,p) odefun(z,p);
        opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
        [~,y0] = ode45(@(t,y) odefw(t,y,p0), tspan, z0, opts);
        [~,y] = ode45(@(t,y) odefw(t,y,p0), [0 2*pi/p0(1)], y0(end,:));
        [~, warnId] = lastwarn;

        if any(isnan(y(:))) || strcmp(warnId,'MATLAB:ode45:IntegrationTolNotMet')
            warning('Reduced dynamics with IRs in polar form diverges with [0.1 0.1 0.1 0.1]');
        else
            z0 = y(end,:)';
        end
end

if ispolar % regularize initial solution if it is in polar form
    z0(2:2:end) = mod(z0(2:2:end),2*pi); % phase angles in [0,2pi]
    for k=1:m
        if z0(2*k-1)<0
            z0(2*k-1) = -z0(2*k-1);      % positive amplitudes
            z0(2*k) = z0(2*k)+pi;
        end
    end
end
