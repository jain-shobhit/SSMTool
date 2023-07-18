% % 

n = 4;
m = 1;
k = 1;
c = 2;
kappa2 = 1;
kappa3 = 0.5;

[M,C,K,fnl,~] = build_model(n,m,c,k,kappa2,kappa3);
%% Dynamical system setup 
DS = DynamicalSystem();
% set(DS,'M',M,'C',C,'K',K,'fnl',fnl); % second order works fine
A = [-K zeros(n);zeros(n) M];
B = [C M;M zeros(n)];
set(DS,'A',A,'B',B,'fnl',{-fnl{1},-fnl{2}});
% set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DS.Options,'Emax',5,'Nmax',100,'notation','tensor')
%% 
epsilon = 5e-3;
f_0 = ones(n,1);
%% 
% Fourier coefficients of Forcing
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();
%% 
S = SSM(DS);
% set(S.Options, 'reltol', 0.1,'notation','multiindex')
set(S.Options, 'reltol', 10,'notation','tensor')
masterModes = [7 8]; % single mode will fail 
order = 5;
S.choose_E(masterModes);
[W_0,R_0] = S.compute_whisker(order);