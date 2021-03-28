%% 
% We extract FRC of a three DOFs model with 1:1:1 internal resonance
% 
% $$\ddot{x}_1+x_1+\epsilon c_1\dot{x}_1+ K(x-y)^3=\epsilon f_1\cos\Omega t,\\\ddot{x}_2+x_2+\epsilon 
% c_2\dot{x}_2+ K[(y-x)^3+(y-z)^3]=\epsilon f_2\cos\Omega t,\\\ddot{x}_3+x_3+\epsilon 
% c_3\dot{x}_3+ K(z-y)^3=\epsilon f_3\cos\Omega t.$$

clear all, close all, clc
%% Example Setup

c1 = 1e-1;
c2 = 2e-1;
c3 = 3e-1;
K = 0.2; 
epsilon = 5e-3;
[mass,damp,stiff,fnl,fext]=build_model(c1,c2,c3,K,epsilon);
%% Dynamical System Setup

DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Forcing
kappas = [1; -1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
%% 
% *Linear Modal analysis*

startMD = tic;
[V,D,W] = DS.linear_spectral_analysis();
timings.MD = toc(startMD);
%% 6D SSM based computation
% *1:1 internal resonance*

S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
resonant_modes = [1 2 3 4 5 6];
mFreq  = [1 1 1];
order = 3;
outdof = [1 2 3];
%% 
% Using ep toolbox - polar (failed when $\rho_3\to 0$

set(S.Options, 'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'method','continuation ep')
set(S.FRCOptions, 'outdof',outdof)
set(S.FRCOptions, 'initialSolver', 'fsolve');
set(S.contOptions,'PtMX',250);

freqRange = [0.995 1.016];
% call extract_FRC to calculate the FRC
FRC = S.extract_FRC('freq',freqRange,order);
%%
% call FRC_cont_ep to calculate the FRC
startep = tic;
FRC_ep_polar = S.FRC_cont_ep('isol_polar',resonant_modes, order, mFreq, 'freq', freqRange);
timings.epPolarFRC = toc(startep);
%% 
% Using ep toolbox - Cartesian

set(S.FRCOptions, 'coordinates','cartesian');
startep = tic;
FRC_ep_cart = S.FRC_cont_ep('isol_cart',resonant_modes, order, mFreq, 'freq', freqRange);
timings.epCartFRC = toc(startep);
%%
% plot results
f1 = figure('Name','Polar-Norm');
f2 = figure('Name',['Polar-Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];
plot_FRC_full(FRC_ep_polar,outdof,order,'freq','lines',figs,'blue')

f3 = figure('Name','Cart-Norm');
f4 = figure('Name',['Cart-Amplitude at DOFs ' num2str(outdof)]);
figs = [f3, f4];
plot_FRC_full(FRC_ep_cart,outdof,order,'freq','lines',figs,'blue')
%% Validation using COCO

nCycles = 100;
coco = cocoWrapper(DS, nCycles, outdof);
set(coco.Options, 'PtMX', 500);
set(coco.Options, 'NAdapt', 1, 'h_max', 50);

startcoco = tic;
bd2 = coco.extract_FRC(freqRange);
timings.cocoFRC = toc(startcoco)