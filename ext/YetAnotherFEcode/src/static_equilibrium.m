function [ u_lin, u ] = static_equilibrium( Assembly, Fext, varargin )
% finds the equilibrium configuration of the model subject to Fext load.
%   Detailed explanation goes here
K = Assembly.DATA.K;
u_lin = Assembly.solve_system(K,Fext);
u0 = Assembly.constrain_vector(u_lin);

[nsteps,tol,method] = parse_inputs(varargin{:});
switch method
    case 'fsolve'
        options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'MaxIterations',10000);
        [ueq] = fsolve(@(u)f(u,Assembly,Fext),u0,options);
        u = Assembly.unconstrain_vector(ueq);
        
    case 'newton'
        u = u_lin/nsteps;
        figure; xlabel('Normalized load');ylabel('$$\|\mathbf{u}\|$$')
        h = animatedline;
        addpoints(h,0,0);
        for j = 1:nsteps
            Fext_j = j*Fext/nsteps;
            c0 = norm(Assembly.constrain_vector(Fext_j));
            it = 0;
            while true
                [K, Fint] = Assembly.tangent_stiffness_and_force(u);
                residual = Fext_j - Fint;
                c = norm(Assembly.constrain_vector(residual))/c0;
                fprintf('STEP %d, ITERATION %d, RESIDUAL %d \n',j,it,c);
                if c < tol
                    break
                end
                correction = Assembly.solve_system(K,residual);
                u = u + correction;
                it = it + 1;
            end
            addpoints(h,j/nsteps,norm(u));
            drawnow 
        end
end

end

function [F,K] = f(u,Assembly,Fext)
x = Assembly.unconstrain_vector(u);
[Kt, Fint] = Assembly.tangent_stiffness_and_force(x);
K = Assembly.constrain_matrix(Kt);
F = Assembly.constrain_vector(Fint - Fext);
end


function [nsteps,tol,method] = parse_inputs(varargin)
%% parsing inputs
defaultnsteps = 100;
defaulttol = 1e-6;
defaultmethod = 'fsolve';

p = inputParser;
addParameter(p,'nsteps',defaultnsteps, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'tol',defaulttol, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'method',defaultmethod,@(x)validateattributes(x, ...
    {'char'},{'nonempty'}))
parse(p,varargin{:});

nsteps = p.Results.nsteps;
tol = p.Results.tol;
method = p.Results.method;
end