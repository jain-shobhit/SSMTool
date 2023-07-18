%% 2 DOF Parametric Duffing Oscillator
% We reproduce results for example of a 2DOF coupled oscillator system with 
% parametric excitation as given in Szabelski, K. & Warminki, J. Vibration of 
% a Non-Linear Self-Excited System with Two Degrees of Freedomunder External and 
% Parametric Excitation.Nonlinear Dynamics 1997 14:114,23–36 (1 1997). <https://doi.org/10.1023/A:1008227315259 
% https://doi.org/10.1023/A:1008227315259>
% 
% 
% 
% The negative linear damping leads to exponentially growing amplitudes, the 
% effect of which is balanced with positive nonlinear damping which starts contributing 
% at significant oscillation amplitudes.

clear all; close all; clc
%% Generate model

[M,C,K,fnl,fext] = build_model();

%% Dynamical system setup
% We consider the forced  and parametrically excited system
% 
% $$\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}+\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})=\epsilon\mathbf{g}(\mathbf{x},\dot{\mathbf{x}},\Omega 
% t),$$
% 
% which can be written in the first-order form as 
% 
% $$\mathbf{B}\dot{\mathbf{z}}	=\mathbf{Az}+\mathbf{F}(\mathbf{z})+\epsilon\mathbf{G}(\mathbf{z},\phi),\\\dot{\mathbf{\phi}}	
% =\mathbf{\Omega}$$
% 
% where
% 
% $$\mathbf{z}=\left[\begin{array}{c}\mathbf{x}\\\dot{\mathbf{x}}\end{array}\right],\quad\mathbf{A}=\left[\begin{array}{cc}-\mathbf{K} 
% & \mathbf{0}\\\mathbf{0} & \mathbf{M}\end{array}\right],\mathbf{B}=\left[\begin{array}{cc}\mathbf{C} 
% & \mathbf{M}\\\mathbf{M} & \mathbf{0}\end{array}\right],\quad\quad\mathbf{F}(\mathbf{z})=\left[\begin{array}{c}\mathbf{-\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})}\\\mathbf{0}\end{array}\right],\quad\mathbf{G}(\mathbf{z},\mathbf{\phi})=\left[\begin{array}{c}\mathbf{g}(\mathbf{x, 
% \dot{x}},\phi)\\\mathbf{0}\end{array}\right]$$

% Dynamical System
DS = DynamicalSystem();
set(DS,'order',2)
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
%% Add forcing
% The dynamical system is forced externally and parametrically with $\mathbf{g(x},\Omega 
% t) = 4 \mu \left[\begin{array}{cc} k & -k \\ -k &  k \end{array}\right] \cos(2\Omega 
% t) \mathbf{x} +\left[\begin{array}{cc}    q \cos (\Omega t )   \\ 0\end{array}\right]$

DS.add_forcing(fext,1e-1);
%% Linear Modal Analysis 

% Analyse spectrum
[V,D,W] = DS.linear_spectral_analysis();
% Choose Master subspace (perform resonance analysis)

S = SSM(DS);
set(S.Options, 'reltol', 0.5,'notation','multiindex')

%Choose Master subspace
masterModes = [3,4];
S.choose_E(masterModes);
%% Forced response curves using SSMs
% Obtaining *forced response curve* in reduced-polar coordinate

order = 5; % Approximation order
%% 
% setup options

outdof = [1,2];
set(S.Options, 'reltol', 0.5,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 80, 'nPar',600, 'nPsi', 80, 'rhoScale', 4 )
set(S.FRCOptions, 'method','level set') % 'continuation ep'
set(S.FRCOptions, 'outdof',outdof)
set(S.FRCOptions,'coordinates','cartesian')
%% 
% choose frequency range around the master mode frequency


omega0 = imag(S.E.spectrum(1));
OmegaRange =omega0*[0.95 1.1];
%% 
% Extract forced response curve

startFRCSSM = tic;
FRC = S.extract_FRC('freq',OmegaRange,order);
figFRC = gcf;
timings.FRCSSM = toc(startFRCSSM);
%% Verification: Collocation using <https://sourceforge.net/p/cocotools/wiki/Home/ coco>
% Dankowicz, H., & Schilder, F. (2013).  _Recipes for Continuation,_ SIAM Philadelphia. 
% <https://doi.org/10.1137/1.9781611972573 <https://doi.org/10.1137/1.9781611972573>>

nCycles = 10;

coco = cocoWrapper(DS, nCycles, outdof);
set(coco,'initialGuess','forward')
set(coco.Options, 'NAdapt', 1);
set(coco.Options,'ItMX',100,'NTST', 70,'PtMX',350); 

figure(figFRC)
hold on
startcoco = tic;
bd3 = coco.extract_FRC(OmegaRange);
timings.cocoFRC = toc(startcoco);