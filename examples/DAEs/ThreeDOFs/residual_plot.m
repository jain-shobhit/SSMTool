% This file reproduces Fig. 6 of the DAE manuscript. This routine should be
% executed after demo_none.

clear all
% Setup model
om1 = 2;
om2 = 3;
om3 = 5;
zeta1 = 0.01;
zeta2 = 0.05;
zeta3 = 0.05;
f1    = 1;
[B,A,Fnl,Fext] = build_model(om1,om2,om3,zeta1,zeta2,zeta3,f1,'none');
DS = DynamicalSystem();
set(DS,'B',B,'A',A,'fnl',Fnl);
set(DS.Options,'Emax',6,'Nmax',20,'notation','multiindex')
epsilon = 0.03;

kappas = [1; -1];
coeffs = [Fext Fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
% 
% *Linear Modal analysis*
[V,D,W] = DS.linear_spectral_analysis();
% Autonomous SSM 
S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
resonant_modes = [1 2]; % choose master spectral subspace
order = 11;                  % SSM expansion order
S.choose_E(resonant_modes);
[W0,R0] = S.compute_whisker(order);
sol = ep_read_solution('none-O11eps3.ep',3); % peak SN
% Non-autonomous SSM
S.System.Omega = sol.p(1);
[W1, R1] = S.compute_perturbed_whisker(order);

nt  = 128;
phi = linspace(0,2*pi,nt);
% reduced coordinate
state = sol.x(1)*exp(1i*sol.x(2));
state = state*exp(1i*phi);
p = [state;conj(state)];
% time history for residual
y = S.compuate_invariance_residual(W0,R0,p,'nonauto',W1,R1);
set(gca,'FontSize',14); grid on;


