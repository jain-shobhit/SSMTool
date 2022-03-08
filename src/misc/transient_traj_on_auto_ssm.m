function traj = transient_traj_on_auto_ssm(DS, modes, W_0, R_0, tf, nsteps, outdof, z0, varargin)
% TRANSIENT_TRAJ_ON_AUTO_SSM This function returns transient trajectory on
% autonomous SSM. Given an initial condition z0 in full system, we first
% project it on SSM (based on linear projection) and then perform time
% integration using ode45 for t\in[0,tf]. Trajectory at outdof is avaliable
% in the output. The inputs W_0 and R_0 can be obtained from compute_whisker
% routine. DS is the dynamical system class and modes represent master
% modes.
%
% varargin = s0 (initial condition in reduced coordinates)

%% project z0 to reduced manifold 
if numel(varargin)==0
    qinit = proj2SSM(z0,'linear',DS.spectrum.W(:,modes),DS.B);
else
    qinit = varargin{1};
end

%% forward simulation in reduced dynamics 
% extract coefficients and exponents
beta  = [];
kappa = [];
for k = 2:numel(R_0)
    R = R_0{k};
    betak = R.coeffs;
    if ~isempty(betak)
        kappak = R.ind;
        % assemble terms
        beta  = [beta betak];
        kappa = [kappa; kappak];
    end
end
autData.lamd  = DS.spectrum.Lambda(modes);
autData.beta  = beta;
autData.kappa = kappa;

% Construct ode45-compatible vector field 
odefun = @(t,z) auto_red_dyn(z,autData);
% forward simulation of reduced dynamics
tsamp  = linspace(0,tf,nsteps+1);
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y]  = ode45(@(t,y) odefun(t,y), tsamp, qinit, options);
% mapping it back to physical domain
state  = transpose(y);
z    = reduced_to_full_traj([],state,W_0);
zout = z(outdof,:);

traj.time = t;
traj.red  = y;
traj.phy  = zout';

end