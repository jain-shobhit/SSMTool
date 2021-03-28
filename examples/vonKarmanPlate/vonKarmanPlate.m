%% A von Karman plate with 1:1 internal resonance
%% Geometry and finite element model

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
%% 4D-SSM based reduction

set(S.FRCOptions, 'omegaSampStyle','cocoBD');
set(S.contOptions, 'PtMX', 2500, 'h_max', 5, 'h_min', 1e-4, 'ItMX', 20);
set(S.FRCOptions, 'initialSolver', 'fsolve', 'outdof', outdof);
resonant_modes = [3 4 5 6];
freqrange = [0.95 1.1]*imag(D(3));
mFreq   = [1 1];
order   = 5;
FRC = S.FRC_cont_ep('isol_polar',resonant_modes, order, mFreq, 'freq', freqrange);
%%
% plot results
f1 = figure('Name','Norm');
f2 = figure('Name',['Amplitude at DOFs ' num2str(outdof')]);
figs = [f1, f2];
plot_FRC_full(FRC,outdof,order,'freq','lines',figs,'blue')