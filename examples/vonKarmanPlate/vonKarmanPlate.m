%% A von Karman plate with 1:1 internal resonance
% In this example, we consider a simply-supported von Karman square plate subject 
% to harmonic excitation. Due to the geometric symmetry, the natural frequenices 
% of the second and the third modes are the same. In other words, the system has 
% 1:1 internal resonance between the two modes. We extract the forced response 
% curve using SSM reduction. 
%% Geometry and finite element model
% We discretize the plate with 200 elements and then 606 degrees-of-freedom. 
% The dimension of the phase space of the full system is then 1212. As we will 
% see, the reduced-order model obtained by SSM reduction is four dimensional.

l = 1; % length of domain [m]
b = 1;  % breadth of domain [m]
t = 1e-2; % thickness of plate [m]
w = 0.0; % curvature parameter (height of the midpoint relative to ends) [m]
% material properties
E     = 70e9;  % 70e9 % 200e9 % Young's modulus [Pa]
rho   = 2700; % 2700 % 7850 % density [kg/m^3]
nu    = 0.33;    % Poisson's ratio 
kappa = 1e5; % material damping modulus 1e8

% Mesh
nElements = 10;
nl = nElements;
nb = nElements;
bc = 'SSSS';
startFE = tic;
[M,C,K,fnl,~,outdof] = build_model(l,b,t,w,E,rho,nu,kappa,bc,nl,nb);
timeFE = toc(startFE);
varargout{1} = timeFE;
n  = length(M);
%% Dynamical System Setup

notation = 'multiindex';
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation',notation);
%% Linear Modal analysis

[V,D,W] = DS.linear_spectral_analysis();
%% Add forcing

f_0 = zeros(n,1);
f_0(outdof(1)) = 100;
kappas  = [-1; 1];
epsilon = 0.5;
coeffs  = epsilon * [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas);
%% Primiary resonance with IRs
% *Create SSM*

S = SSM(DS);
set(S.Options, 'reltol', 0.8,'notation',notation); % 0.8 is used
%% 
% *Extract FRC*
% 
% Recall that the natural frequencies of the linear modes are given by
% 
% $$\omega_{(i,j)}=\left(\frac{i^2}{l^2}+\frac{j^2}{b^2}\right)\pi^2\sqrt{\frac{D}{\rho 
% t}}$$
% 
% where $D$ is the bending stiffness. Given $l=b$, we have $\omega_{(1,2)}=\omega_{(2,1)}$, 
% namely the 1:1 internal resonance between the second and the third bending modes. 
% In addition, we have $\omega_{(3,1)}=\omega_{(1,3)}=2\omega_{(1,2)}$, which 
% implies that the higher modes (1,3) and (3,1) will also be taken as a part of 
% the spectral subspace in the extract_FRC routine. Here we are mainly concerned 
% with the 1:1 internal resonance. So we need to specify the spectral subspace, 
% which can be easily done in the SSM-ep toolbox.

set(S.FRCOptions, 'omegaSampStyle','cocoBD','method','continuation ep');
set(S.contOptions, 'PtMX', 2500, 'h_max', 5, 'h_min', 1e-4, 'ItMX', 20);
set(S.FRCOptions, 'initialSolver', 'fsolve');
freqrange = [0.95 1.1]*imag(D(3));
order     = 5;
resonant_modes = [3 4 5 6];
mFreq   = [1 1];
FRC = S.SSM_isol2ep('isol_polar',resonant_modes, order, mFreq, 'freq', freqrange,outdof);
%% 
% *Exercises*
%% 
% # Use SSM_isol2ep to calculate the forced response curve with varied forcing 
% amplitude $\epsilon$.
% # Use SSM_ep2SN to locate the solution manifold of saddle-node bifurcation 
% periodic orbits.
% # Use SSM_ep2HB to locate the solution manifold of torus bifurcation periodic 
% orbits.