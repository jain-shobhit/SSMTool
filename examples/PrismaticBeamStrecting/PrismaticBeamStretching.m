%% 
% We consider a clamped-pinned beam. With such a boundary condition, the first 
% two modes have near 1:3 internal resonance as follows
% 
% $$\omega_2=3\omega_1(1+\epsilon\sigma_1)$$
% 
% with $\epsilon\sigma_1=0.0801$.
% 
% Nayfeh [1] inverstigated the forcec response of such a system under external 
% harmonic response. Specifically, modal expansion (with linear modes) is used 
% to transfer PDEs to a set of ODEs
% 
% $$\ddot{u}_n+\omega_n^2u_n=\epsilon\left(-2c_n\dot{u}_n+F_n\cos\lambda t+\nu\sum_{m,p,q}\alpha_{nmpq}u_mu_pu_q\right), 
% n=1,\cdots,$$
% 
% Then multiple scale perturbation method is used to study the forced response 
% curve. Interestingly, when the excitation frequency is around the second natural 
% frequency, there exits solution branches where the amplitude of the first mode 
% dominiates the system response. Here we use SSM reduction to study such a system.
% 
% [1] Nayfeh, A. H., Mook, D. T., & Sridhar, S. (1974). Nonlinear analysis of 
% the forced response of structural elements. _The Journal of the Acoustical Society 
% of America_, _55_(2), 281-291.
%% Setup Dynamical System

clear all;
epsilon = 1e-4;
c  = 100;
f = 5/epsilon;
n = 10;               % number of modes
Fext = zeros(n,1);    % excitation at modal coordinate
Fext(1) = f;
[mass,damp,stiff,fnl,fext] = build_model(c,Fext,epsilon,n);

% Create 
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% Forcing
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas);
%% Linear Modal analysis

[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
resonant_modes = [1 2 3 4];
order = 3;
outdof = [1 2];
%% Primary resonance of the first mode

freqrange = [0.98 1.06]*imag(D(1));
set(S.FRCOptions, 'nCycle',500, 'initialSolver', 'forward');
set(S.contOptions, 'PtMX', 300, 'h_max', 0.5);
set(S.FRCOptions,'outdof', outdof, 'omegaSampStyle', 'cocoBD');
start = tic;
FRC_ep1 = S.FRC_cont_ep('isol',resonant_modes, order, [1 3], 'freq', freqrange);
timings.FRC_ep1 = toc(start)
% Load Nayfeh's solution for comparison

f1 = figure('Name','P1-Norm');
f2 = figure('Name',['P1-Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];
plot_FRC_full(FRC_ep1,outdof,order,'freq','lines',figs,'blue')


Nayfeh1 = load('NayfehFirstMode');
om1 = 3.8553;
omsamp = om1*(1+Nayfeh1.epsilon*Nayfeh1.SIG2);
figure(f2); hold on
subplot(2,1,1)
plot(omsamp(1:100:end),Nayfeh1.A1(1:100:end),'ro','DisplayName','Nayfeh');
subplot(2,1,2)
plot(omsamp(1:100:end),Nayfeh1.A2(1:100:end),'ro','DisplayName','Nayfeh');

% Nonzero forcing at other modes $f_1=f_2=\cdots=f_{10}\neq0$

Fext = zeros(n,1)+f;  % excitation at modal coordinate
[~,~,~,~,fext2] = build_model(c,Fext,epsilon,n);
coeffs = [fext2 fext2]/2;
DS.add_forcing(coeffs, kappas);
outdof2 = [2 3]

set(S.FRCOptions,'outdof',outdof2)
start = tic;
FRC_ep2 = S.FRC_cont_ep('isolAllf',resonant_modes, order, [1 3], 'freq', freqrange);
timings.FRC_ep2 = toc(start)
%%
f1 = figure('Name','P1allF-Norm');
f2 = figure('Name',['P1allF-Amplitude at DOFs ' num2str(outdof2)]);
figs = [f1, f2];
plot_FRC_full(FRC_ep2,outdof2,order,'freq','lines',figs,'blue')
 
%% Primary resonance of the second mode
% _*epsilon = 1e-4; c = 10; f2 = 40/epsilon;*_ 

c  = 10;
f = 40/epsilon;
n = 10;               % number of modes
Fext = zeros(n,1);    % excitation at modal coordinate
Fext(2) = f;
[mass,damp,stiff,fnl,fext] = build_model(c,Fext,epsilon,n);

% Create dynamical system
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% Forcing
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas);
%% Linear Modal analysis

[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 1.0,'notation','multiindex');
resonant_modes = [1 2 3 4];
order = 3;
outdof = [1 2];
freqrange = [0.94 1.12]*imag(D(3));

set(S.FRCOptions, 'omegaSampStyle', 'cocoBD');
set(S.contOptions, 'h_max', 5, 'PtMX', 250);
set(S.FRCOptions, 'initialSolver', 'fsolve');
set(S.FRCOptions, 'coordinates', 'cartesian'); 
%% 
% *Non-zero* $a_1$ *(reproduce Nayfehs' solutions to get initial guess)*

p0 = [imag(D(3));1];
z0 = [97.8562  -55.4066    1.9095    0.9820]'*1e3;
set(S.FRCOptions,'p0',p0,'z0',z0);
set(S.FRCOptions,'outdof',outdof);
start = tic;
FRC_ep21 = S.FRC_cont_ep('isol2',resonant_modes, order, [1/3 1], 'freq', freqrange);
timings.FRC_ep21 = toc(start)
% Load Nayfeh's solution for comparison

f1 = figure('Name','P2nonzero-Norm');
f2 = figure('Name',['P2nonzero-Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];
plot_FRC_full(FRC_ep21,outdof,order,'freq','lines',figs,'blue')


Nayfeh21 = load('NayfehSecondModeA1Nequal0.mat');
om2 = 12.4927;
omsamp = om2*(1+1e-4*Nayfeh21.SIG2);
figure(gcf); hold on
subplot(2,1,1)
plot(omsamp(1:150:end),Nayfeh21.A1(1:150:end),'ro','DisplayName','Nayfeh');
subplot(2,1,2)
plot(omsamp(1:150:end),Nayfeh21.A2(1:150:end),'ro','DisplayName','Nayfeh');
%% 
% *Zero* $a_1$ *upper branch*

set(S.FRCOptions,'p0',[],'z0',[]);
start = tic;
FRC_ep22 = S.FRC_cont_ep('isol3',resonant_modes, order, [1/3 1], 'freq', freqrange);
timings.FRC_ep22 = toc(start)
%% 
% *Zero* $a_1$ *lower branch*

p0 = [1.06*imag(D(3));1];
z0 = [0 0 1 1]';
set(S.FRCOptions,'p0',p0,'z0',z0);
start = tic;
FRC_ep23 = S.FRC_cont_ep('isol4',resonant_modes, order, [1/3 1], 'freq',freqrange);
timings.FRC_ep23 = toc(start)
% Plot zero $a_1$ results in the same figure

f1 = figure('Name','P2zero-Norm');
f2 = figure('Name',['P2zero-Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];
plot_FRC_full(FRC_ep22,outdof,order,'freq','lines',figs,'blue')
plot_FRC_full(FRC_ep23,outdof,order,'freq','lines',figs,'blue')

% Load Nayfeh's solution for comparison

Nayfeh21 = load('NayfehSecondModeA1equal0.mat');
om2 = 12.4927;
omsamp = om2*(1+Nayfeh21.epsilon*Nayfeh21.Sig2);
figure(gcf); hold on
subplot(2,1,2)
plot(omsamp(1:100:end),Nayfeh21.A2(1:100:end),'ro','DisplayName','Nayfeh');