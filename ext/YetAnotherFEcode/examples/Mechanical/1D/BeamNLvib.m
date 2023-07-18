%% EXAMPLE 1D-BEAM using NLvib
% Example based on:
% Sombroek et al. (2018). Numerical computation of nonlinear normal modes 
% in a modal derivative subspace. Computers and Structures, 195, 34–46. 
% https://doi.org/10.1016/j.compstruc.2017.08.016
%
% Author: Jacopo Marconi, Politecnico di Milano
% Created: 21 April 2021
% Last modified: 27 April 2021

clear
clc
close all

imod = 1; % mode analyze

% NOTE: you can load "BeamNLvib.mat" which contains the results for the
% beam meshed with 8 elements (all other parameters set as in this script).
% Run the PLOT sections to inspect the results.

%% Parameters                                                       
% geometry
len = 1;        	% length
height = 1e-2;    	% height in the bending direction
thickness = 1e-2;	% thickness in the third dimension

% mesh
nElements = 4;
dx = len / nElements;

% material properties
E       = 210e9;    % Young's modulus
rho     = 7800;     % density
nu      = 0.3;      % nu

%% Structural Model                                                 
% Material
myBeamMaterial  = KirchoffMaterial();
set(myBeamMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);

% Element
myElementConstructor = @()BeamElement(thickness, height, myBeamMaterial); % same element all across the domain

% Mesh
x = (0:dx:len).';
nNodes = size(x,1);
Nodes = [x, zeros(nNodes,1)];

Elements = [1:nNodes-1;2:nNodes].';
BeamMesh = Mesh(Nodes);
BeamMesh.create_elements_table(Elements,myElementConstructor);

nNodes = BeamMesh.nNodes;
BeamMesh.set_essential_boundary_condition([1 nNodes],1:3,0);

ndofs = length( BeamMesh.EBC.unconstrainedDOFs );
ntot = BeamMesh.nDOFs;

% ASSEMBLY ________________________________________________________________
BeamAssembly = Assembly(BeamMesh);
M = BeamAssembly.mass_matrix();
u0 = zeros( BeamMesh.nDOFs, 1);
[K,~] = BeamAssembly.tangent_stiffness_and_force(u0);

% store matrices
BeamAssembly.DATA.K = K;
BeamAssembly.DATA.M = M;

% External force __________________________________________________________
forced_dof = round(nNodes/2) * 3 - 1; % middle node, vertical dof
Fext = zeros(ntot, 1);
Fext( forced_dof ) = 1;
Fextc = BeamAssembly.constrain_vector( Fext );

% Let us also define the index of the forced dof in the constrained vector:
forced_dof_c = BeamAssembly.free2constrained_index( forced_dof );

%% Eigenvalue problem                                               
n_VMs = ndofs; % first n_VMs modes with lowest frequency calculated 
Kc = BeamAssembly.constrain_matrix(K);
Mc = BeamAssembly.constrain_matrix(M);
[Phi,om2] = eigs(Kc, Mc, n_VMs, 'SM');
[om, ind] = sort(sqrt(diag(om2)));
f0 = om/2/pi;
Phi = Phi(:,ind);
for ii = 1:n_VMs
    Phi(:,ii) = Phi(:,ii)/max(sqrt(Phi(1:3:end,ii).^2+Phi(2:3:end,ii).^2));
end
Phi = BeamAssembly.unconstrain_vector(Phi);

u = Phi(1:3:end, imod);
v = Phi(2:3:end, imod);
x = Nodes(:, 1);
y = Nodes(:, 2);

figure
scale = len / 3;
plot(x, y, 'k--o'); hold on
plot(x+u*scale, y+v*scale, 'b-o')
grid on; axis equal tight; 
xlabel('x [m]'); ylabel('y [m]'); title(['mode ' num2str(imod)])
drawnow

% Damping _________________________________________________________________
Qfactor = 100;
csi = 1./(2*Qfactor);   % dimensionless damping
om0 = 2*pi*f0(1);       % first eigenfrequency
alfa = 2 * om0 * csi;
D = alfa*M;             % Mass-proportinal damping: D = a*M
BeamAssembly.DATA.D = D;
Dc = BeamAssembly.constrain_matrix(D);



%% (1) NLvib: NMA - Harmonic Balance                                
% Example adapted from "09_beam_cubicSpring_NM" from NLvib

BeamSystem = FE_system( BeamAssembly, Fext );

PHI_lin = BeamAssembly.constrain_vector(Phi);

% Analysis parameters _____________________________________________________
H = 9;              % harmonic order (size of the problem is (2H+1)*ndof,
                    % 	    i.e. A0 and (Ai,Bi) for i=1...H for each dof)
N = 3*H+1;       	% number of time samples per period
log10a_s1 = -5;     % start vibration level (log10 of modal mass)
log10a_e1 = -2;     % end vibration level (log10 of modal mass)
inorm = ndofs-1;   	% coordinate for phase normalization

% INITIALIZATION __________________________________________________________
% Initial guess vector x0 = [Psi; om; del], where del is the modal
% damping ratio, estimate from underlying linear system
omi = om(imod);                 % linear eigenfrequency
phi = PHI_lin(:,imod);          % linear eigenmode
Psi = zeros((2*H+1)*ndofs, 1);  % initialize ALL the harmonics to zero
Psi(ndofs+(1:ndofs)) = phi;     % initialize the first harmonic with the 
                                % linear eigenmode. Psi(1:ndofs) are the
                                % 0th harmonic terms (static).
x0 = [Psi; omi; 0];          	% initial guess
q_scl = max(abs(Psi));        	% used for scaling (scl)

% Solve and continue w.r.t. Om ____________________________________________
ds = .02;                       % arclength parameter
Sopt=struct('dynamicDscale',1);	% set scaling and other options here
f_scl = mean(abs(Kc*phi));      % scaling of dynamic force equilibrium
fun_postprocess = {};           % (optional) feval(fun_postprocess, X),
                                %            results are stored in Sol
% solve
[X1, Solinfo, Sol] = solve_and_continue(x0,...
    @(X) HB_residual(X, BeamSystem, H, N, 'NMA', inorm, f_scl),...
    log10a_s1, log10a_e1, ds, [Sopt, fun_postprocess]);

% Interpret solver output
r1 = nlvib_decode(X1, Solinfo, Sol, 'NMA', 'HB', ndofs, H);
results.NMA.HB = r1;

% POSTPROCESSING __________________________________________________________
% [adapted from EXAMPLE 9 of NLvib]
% Determine total energy in the system from the displacement and velocity
% at t=0

energy = zeros(size(r1.a));
k = repmat(1:H, ndofs, 1); % harmonic numbers

for ii = 1 : length( r1.omega )
    % q(t) = real(sum_k(Qc*exp(i*k*W*t))), q0 = q(t=0)
    % u(t) = dq/dt = real(sum_k(Qc*(i*k*W)*exp(i*k*W*t))), u0 = u(t=0)
    q0 = r1.Q0(:, ii) + sum(squeeze(r1.Qre(:, ii, :)), 2);
    u0 = sum(squeeze(r1.Qim(:, ii, :)) .* k * r1.omega(ii), 2);
    
    % the total energy is constant on the orbit (therefore we compute it
    % for t=0):
    energy(ii) = 1/2*u0'*Mc*u0 + 1/2*q0'*Kc*q0;
    % ... we should add +1/3*K3*q0^3 +1/4*K4*q0^4, but we don't have K3,K4
end
fprintf(['\nBE WARE: the energy computed here is not correct, only \n', ...
    'convenient (we don''t have third and fourth order tensors K3 \n', ...
    'and K4 of the full model to compute 1/3*K3*q^3 and 1/4*K4*q^4).\n\n'])

results.NMA.HB.energy = energy;

%% (1) PLOT                                                         

r1 = results.NMA.HB;
wp = r1.omega;

figure('units', 'normalized', 'position', [.33 .1 .33 .8])
subplot 311
semilogy(wp, r1.a, 'k-');       % Amplitude vs frequency
    xlabel('\omega [rad/s]'); ylabel('a_{HB}');
    grid on; box on; axis tight
    hold on
    plot(xlim, 10^log10a_s1*[1 1],'r--')
    plot(xlim, 10^log10a_e1*[1 1],'r--')
    ylim(10.^[log10a_s1-.5 log10a_e1+.5])
    title('NMA with Harmonic Balance')
subplot 312
semilogx(r1.energy, wp, 'b-')	% Energy vs frequency (linear)
    grid on; axis tight
    xlabel('energy [J]'); ylabel('\omega [rad/s]');
    ylim([320 500])
    xlim([1e-4 1e1])
    title('Energy-Frequency plot')
subplot 313
h = 1;
A = r1.Qre(forced_dof_c, :, h);
B = r1.Qim(forced_dof_c, :, h);
a1 = squeeze(sqrt( A.^2 + B.^2 ));
plot(wp, a1 / height, 'linewidth', 1);
    xlim([1e-4 1e1])
    xlabel('\omega [rad/s]'); ylabel('|Q_1| / height [-]');
    grid on; box on; axis tight
    title('1^{st} harmonic')
drawnow



%% (2) NLvib:  FRF - Harmonic Balance                             	

omi = om(imod);         % linear eigenfrequency

% Analysis parameters
H = 7;                  % harmonic order
N = 3*H+1;              % number of time samples per period
Om_s = omi * 0.95;   	% start frequency
Om_e = omi * 1.1;    	% end frequency
ds = 1;                 % Path continuation step size
exc_lev = [.5 1];       

fprintf('\n\n FRF from %.2f to %.2f rad/s \n\n', Om_s, Om_e)

% Excitation levels
r2 = cell(length(exc_lev),1);
for iex = 1 : length(exc_lev)
    % Set excitation level
    BeamSystem.Fex1 = Fextc * exc_lev(iex);
    
    % Initial guess (solution of underlying linear system)
    Q1 = (-Om_s^2*Mc + 1i*Om_s*Dc + Kc) \ Fextc;
    y0 = zeros( (2*H+1)*ndofs , 1);
    y0( ndofs + (1:2*ndofs) ) = [real(Q1);-imag(Q1)];
    
    % stuff for scaling
    qscl = max(abs((-omi^2*Mc + 1i*omi*Dc + Kc) \ Fextc));
    dscale = [y0*0+qscl; omi];
    Sopt = struct('Dscale', dscale, 'dynamicDscale', 1);
    
    % Solve and continue w.r.t. Om	
    [X2, Solinfo, Sol] = solve_and_continue(y0, ...
        @(X) HB_residual(X, BeamSystem, H, N, 'FRF'), ...
        Om_s, Om_e, ds, Sopt);
    
    % Interpret solver output
    r2{iex} = nlvib_decode(X2, Solinfo, Sol, 'FRF', 'HB', ndofs, H);
    
    results.FRF.HB{iex} = r2{iex};
end

%% (2) PLOT                                                         

r2 = results.FRF.HB;

figure
h = 1;
for iex = 1 : length(exc_lev)
    % 1st harmonic amplitude of the forced dof (use force_dof_c!)
    A = r2{iex}.Qre(forced_dof_c, :, h);
    B = r2{iex}.Qim(forced_dof_c, :, h);
    W = r2{iex}.omega;
    a1 = squeeze(sqrt( A.^2 + B.^2 ));
    plot(W, a1 / height, 'linewidth', 1); hold on
end
grid on
axis tight
xlabel('\omega [rad/s]')
ylabel('|Q_1| / height [-]')
title('FRF with Harmonic Balance')

try
    % LINEAR RESPONSE
    % compute the linear FRF for comparison
    nw = 201;
    w_linear = linspace(Om_s, Om_e, nw);
    for iex = 1 : length(exc_lev)
        fr_linear = zeros(nw, ndofs);
        for ii = 1:nw
            w = w_linear(ii);
            fr_linear(ii,:) = (-w^2*Mc + 1i*w*Dc + Kc) \ Fextc * exc_lev(iex);
        end
        plot(w_linear, abs(fr_linear(:, forced_dof_c))/height, 'k--')
    end
    drawnow
end


%% (3) NLvib: NMA - shooting                                        

inorm = ndofs-1;                % coordinate for phase normalization
omi = om(imod);                 % linear eigenfrequency
phi = PHI_lin(:,imod);       	% mode shape imod
nnorm = setdiff(1:ndofs,inorm);	% number of remaining modes
qs = phi(nnorm)/phi(inorm);   	% q starting guess (remaining modes only)
                                % (q = displacements)
us = 0*phi(nnorm);            	% u starting guess (remaining modes only)
                                % (u = velocities)
x0 = [qs; us; omi; 0];         	% initial guess

log10a_s2 = -5;   	% start vibration level (as log10 of qs(inorm))
log10a_e2 = -3;    	% end   vibration level (as log10 of qs(inorm))

% NOTE: q(imod) = a = 10.^X(end,:)
%       u(imod) = 0
qscl = max([10^log10a_s2, 10^log10a_e2]);	% scaling factor

% Solve and continue w.r.t. Om ____________________________________________
Np = 1;                         % seek period-Np solutions
Ns = 2^08-1;                   	% time samples for integration
ds = 1e-2;                    	% Path continuation step size 
Sopt = struct('stepmax', 1e4, 'dsmax', ds*10, 'dsmin', ds/1000, ...
            'dynamicDscale', 1);

% *** PostProcessing function
% from X we get the INITIAL CONDITIONS (ICs) in displacement and speed for
% the orbit. If we want to compute the harmonics, we need to run a 
% simulation over the period using these ICs. We can the perform an fft to
% compute the spectrum. This is done in the following function, which also
% returns the DOF-histories in time (displ and speed):
H = 7;
fun_postprocess = @(X) nlvib_ShootingPostprocess(X, BeamSystem, ...
    3*H+1, Np, 'NMA', H, inorm); % this function is run at each solution point

[X3, Solinfo, Sol] = solve_and_continue(x0, ...
    @(X) shooting_residual(X, BeamSystem, Ns, Np, 'NMA', qscl, [], inorm), ...
    log10a_s2, log10a_e2, ds, Sopt, fun_postprocess);

% Interpret solver output _________________________________________________
r3 = nlvib_decode(X3, Solinfo, Sol, 'NMA', 'SH', ndofs, H);
results.NMA.SH = r3;

%% (3) PLOT                                                         

r3 = results.NMA.SH;

fig3 = figure('units', 'normalized', 'position', [.33 .1 .33 .8]);
subplot 311
semilogy(r3.omega, r3.a, '.-')
    xlabel('\omega [rad/s]')
    ylabel('a_{sh}')
    grid on; hold on;
    plot(xlim, 10^log10a_s2*[1 1],'r--')
    plot(xlim, 10^log10a_e2*[1 1],'r--')
    ylim(10.^[log10a_s2-.5 log10a_e2+.5])
    title('NMA with Shooting method')
subplot 312
plot(r3.omega, r3.harmonics.rms.displ(forced_dof_c, :) / height)
    xlabel('\omega [rad/s]')
    ylabel('RMS / height [-]')
    grid on
subplot 313
plot(r3.omega, r3.harmonics.Qabs.displ(forced_dof_c, :, 1+1) / height, 'b')
    xlabel('\omega [rad/s]')
    ylabel('|Q_1| / height [-]')
    grid on
drawnow



%% (4) NLvib: FRF - shooting                                     	

omi = om(imod);   	% mode frequency imod

% Analysis parameters _____________________________________________________
Om_s =  .90* omi; 	% start frequency
Om_e = 1.05* omi;   	% end frequency
Ns = 2^8-1;       	% Number of time samples per period
fprintf('\n\n FRF from %.2f to %.2f rad/s \n\n', Om_s, Om_e)

BeamSystem.Fex1 = Fextc * exc_lev(end);

% Initial guess (solution of underlying linear system) ____________________
Q1 = (-Om_s^2*Mc + 1i*Om_s*Dc + Kc)\BeamSystem.Fex1;
y0 = [real(Q1); -Om_s*imag(Q1)];

% Solve and continue w.r.t. Om ____________________________________________
ds = .1;
Sopt = struct('Dscale',[1e0*ones(size(y0)); Om_s],'dynamicDscale',1);
Np = 1; % we seek period-one solutions

% postprocessing function: evaluates solution in time, computes spectra (up
% to H-th harmonic) and performs stability analysis.
H = 7;
fun_postprocess = @(X) nlvib_ShootingPostprocess(X, BeamSystem, ...
    Ns, Np, 'FRF', H);

[X4, Solinfo, Sol] = solve_and_continue(y0,...
    @(X) shooting_residual(X, BeamSystem, Ns, Np, 'FRF'),...
    Om_s, Om_e, ds, Sopt, fun_postprocess);

r4 = nlvib_decode(X4, Solinfo, Sol, 'FRF', 'SH', ndofs, H);
results.FRF.SH = r4;

%% (4) PLOT                                                         

r4 = results.FRF.SH;

wp = r4.omega;
Q1 = r4.harmonics.Qabs.displ(:, :, 1+1);   % first harmonic
q1 = squeeze(Q1(forced_dof_c, :));
stable = boolean(r4.floquet.stable);

if exist('fig3','var') && isvalid(fig3)
    figure(fig3)
    subplot 313; hold on
else
    figure;
end
q1u = q1; q1u( stable) = NaN;
plot(wp, q1 / height, 'b.-', 'linewidth', 1)
hold on;
plot(wp, q1u / height, 'r.-', 'linewidth', 1)
grid on
xlabel('\omega [rad/s]')
ylabel('Q_1 / height [-]')
title('FRF with Shooting method')


