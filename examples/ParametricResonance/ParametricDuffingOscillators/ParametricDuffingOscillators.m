%% Parametric Duffing Oscillator
% Coupled parametric duffing oscillators along the lines of 
% 
% Szabelski, K. & Warminski, J. Self-excited System Vibrations with Parametric 
% and External Excitations.Journalof Sound and Vibration187,595–607 (4 1995). 
% <https://doi.org/10.1006/jsvi.1995.0547 https://doi.org/10.1006/jsvi.1995.0547> 
% 
% The model considered here is given as 
% 
% 
% 
% The linear and nonlinear stiffness of the coupling spring is varying in a 
% periodic fashion. Furthermore the damping coefficient $\tilde{c}$ is negative, 
% which leads to self excitation. The resulting exponentially growing amplitudes 
% are stabilised by nonlinear damping.

clear all;clc
%% system parameters

mu1 = 0.3;
mu3 = 0.3;
epsilon = 0.2;
%% Generate model

[ M,C,K,fnl,fext] = build_model(mu1,mu3);
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
% $\mathbf{z}=\left[\begin{array}{c}\mathbf{x}\\\dot{\mathbf{x}}\end{array}\right],\quad\mathbf{A}=\left[\begin{array}{cc}-\mathbf{K} 
% & \mathbf{0}\\\mathbf{0} & \mathbf{M}\end{array}\right],\mathbf{B}=\left[\begin{array}{cc}\mathbf{C} 
% & \mathbf{M}\\\mathbf{M} & \mathbf{0}\end{array}\right],\quad\quad\mathbf{F}(\mathbf{z})=\left[\begin{array}{c}\mathbf{-\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})}\\\mathbf{0}\end{array}\right],\quad\mathbf{G}(\mathbf{z},\mathbf{\phi})=\left[\begin{array}{c}\mathbf{g}(\mathbf{x, 
% \dot{x}},\phi)\\\mathbf{0}\end{array}\right]$.

DS = DynamicalSystem();
set(DS,'M',M, 'C', C, 'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
%% Add forcing
% The dynamical system is forced externally and parametrically with $\mathbf{g(x},\Omega 
% t) = \mu\left[\begin{array}{cc} k & -k \\ -k &  k \end{array}\right] \cos(2\Omega 
% t) \mathbf{x} + \kappa\mu \cos (2 \Omega t) \left[\begin{array}{cc}  -(x_1-x_2)^3 
% \\   (x_1-x_2)^3   \end{array}\right]     + \left[\begin{array}{cc}   0 \\ q 
% \cos (\Omega t)   \end{array}\right]$

DS.add_forcing(fext,epsilon);
%% Linear Modal Analysis 

% Analyse spectrum
[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.5,'notation','multiindex')

%Choose Master subspace
masterModes = [1,2];
S.choose_E(masterModes);
%% Forced response curves using SSMs
% Obtaining *forced response curve* in reduced-polar coordinate

order = 5; % Approximation order
%% 
% setup options

outdof = [1:4];
set(S.Options, 'reltol', 0.5,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 10, 'nPsi', 200, 'rhoScale', 6 )
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