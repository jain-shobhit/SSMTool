%% 
% We compute SSM for a Mathieu Equation with nonlinearities
% 
% % 
clear all; close all; clc
%run ../../install.m
% change to example directory
%exampleDir = fileparts(mfilename('fullpath'));
%cd(exampleDir)

%% Parameters
c = 0.1;  % Damping
w0 = 1; % Natural Frequency
lam3 = -0.1;
gam3 = -0.1;
eps  = 0.2;
dampnl = 0;
[Fext,Fnl] = build_model(w0^2,lam3,gam3,dampnl);


%% System of type 
%  B d(x;dot(x))/dt = A(x;dot(x)) +Fnl(z) + Fext(t,z)
B = eye(2);
A = [0,1;-(w0^2),-c];

% Dynamical System
DS = DynamicalSystem();
set(DS,'order',1)
set(DS,'A',A,'B',B,'F',Fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Analyse spectrum
[V,D,W_evec] = DS.linear_spectral_analysis();

% External forcing
DS.add_forcing(Fext,eps);

% Set up SSM object
S = SSM(DS);
set(S.Options, 'reltol', 0.5,'notation','multiindex')

%Choose Master subspace
resModes = [1,2];
S.choose_E(resModes);

% Choose expansion order
order =12;

outdof = 1;


%% Set options
set(S.Options,    'reltol', 0.5,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 80, 'nPar', 80, 'nPsi', 80, 'rhoScale', 4 )
set(S.FRCOptions, 'method','continuation po') % 'continuation ep'
set(S.FRCOptions, 'outdof',outdof)
set(S.FRCOptions, 'coordinates','cartesian')
set(S.FRCOptions, 'resType','2:1') %Forcing at twice the natural frequency 
set(S.FRCOptions, 'periodsRatio',2) %seeking 2:1 periodic response wrt forcing
set(S.FRCOptions, 'branchSwitch','false') %continue BPs of primary branch
set(S.FRCOptions, 'branchSwitchbidir','true') %continue BPs of primary branch

%set(S.contOptions,'h_max',0.05) %continuation step size

%% choose frequency range
omega0 = imag(S.E.spectrum(1));

OmegaRange =[1.85,2.17]; %Subharmonic resonance

%% extract forced response curve

startFRCSSM = tic;
FRC = S.extract_FRC('freq',OmegaRange,order);
figFRC = gcf;
timings.FRCSSM = toc(startFRCSSM);



%% Get results from full system
% {
nCycles = 10;

coco = cocoWrapper(DS, nCycles, outdof);
set(coco,'initialGuess','forward')
set(coco,'branchSwitch','false','periodsRatio',2) % include new branches, 2T periodic response
set(coco.Options, 'NAdapt', 1);
set(coco.Options,'ItMX',100,'NTST', 70,'PtMX',30); %for convergence, smaller stepsize

figure(figFRC)
hold on
startcoco = tic;
bd3 = coco.extract_FRC(OmegaRange);
timings.cocoFRC = toc(startcoco);
%}