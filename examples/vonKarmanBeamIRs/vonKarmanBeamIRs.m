%% In this example, we consider a von Karman beam with linear support spring
% at its midspan. The stiffness of the support spring is tund such that 1:3
% internal resonance occurs.

clear all
%% Model Setup
nElements = 100;
[M,~,K,fnl,~,Outdof] = build_model(nElements);
outdof = Outdof(2); % the point at mid span
n = length(M);
Kc = K;
%% 
% _*1:3 near resonance is observed at*_ $k=37$
%% Dynamical System Setup
kLinear = 37;
Kc(outdof,outdof) = K(outdof,outdof)+kLinear;   % linear part
kNonlinear = 0;                                 % nonlinear part
f3_new  = fnl{2};
fnl_new = fnl;
f3_new(outdof,outdof,outdof,outdof) = f3_new(outdof,outdof,outdof,outdof)+kNonlinear; % cubic nonlinerity
fnl_new{2} = f3_new;
C = (2e-4)/9*Kc;

DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',Kc,'fnl',fnl_new);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');

%% Primiary resonance with IRs
% *Add forcing*
f_0 = zeros(n,1);
f_0(outdof) = 1000;
kappas = [-1; 1];
epsilon   = 0.02;
coeffs = epsilon * [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas);
%% 
% *Create SSM*
S = SSM(DS);
set(S.Options, 'reltol', 0.8,'notation','multiindex')
%% 
% *Using ep toolbox*
om = eigs(Kc,M,2,'smallestabs');
om = sqrt(om);
freqRange = [0.96 1.05]*om(1);
order = 7;
set(S.FRCOptions, 'nCycle',50000);
set(S.FRCOptions, 'initialSolver','forward');
set(S.contOptions, 'h_max', 0.2, 'PtMX', 1000, 'ItMX',20);
set(S.FRCOptions, 'omegaSampStyle','cocoBD');
set(S.FRCOptions, 'method','continuation ep','outdof',Outdof');
FRC = S.extract_FRC('freq', freqRange, order);
