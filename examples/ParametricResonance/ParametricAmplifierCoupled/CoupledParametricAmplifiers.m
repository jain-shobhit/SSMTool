%% Parametric Duffing Oscillator
% The model features two oscillators, each of which resembels the model of a 
% parametric amplifier as presented by 
% 
% Li, D. & Shaw, S. W. The effects of nonlinear damping on degenerate parametric 
% amplification.NonlinearDynamics 2020 102:4102,2433–2452 (4 Dec. 2020) <https://doi.org/10.1006/jsvi.1995.0547 
% https://doi.org/10.1006/jsvi.1995.0547> 
% 
% The model is extended to a two- dimensional setting. 
% 
% 
% 
% The parametric excitation of the linear stiffness is of the form $k(t) = k(1 
% + \mu \cos (\Omega t))$.
%% system parameters

mu = 0.12;
psi = pi/4;
%% Generate model

[ M,C,K,fnl,fext] = build_model(psi,mu);
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

% Dynamical System
DS = DynamicalSystem();
set(DS,'M',M, 'C', C, 'K',K,'fnl',fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')


%% Add forcing
% The dynamical system is forced externally and parametrically with $\mathbf{g(x},\Omega 
% t) = \mu\left[\begin{array}{cc} k & 0 \\ 0 &  k \end{array}\right] \cos(2\Omega 
% t) \mathbf{x} +\left[\begin{array}{cc}   0 \\ q \cos (\Omega t + \psi)   \end{array}\right]$
% 
% Parameters are chosen as $\epsilon = 1, q = 0.2$

% External forcing
DS.add_forcing(fext,1);
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

outdof = [1,2];
set(S.Options, 'reltol', 0.5,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 200, 'rhoScale', 6 )
set(S.FRCOptions, 'method','level set') % 'continuation ep'
set(S.FRCOptions, 'outdof',outdof)

%% 
% choose frequency range around the master mode frequency

omega0 = imag(S.E.spectrum(1));
OmegaRange =omega0*[0.9 1.2];
%% 
% Extract forced response curve

startFRCSSM = tic;
FRC  = S.extract_FRC('freq',OmegaRange,order);
timings.FRCSSM = toc(startFRCSSM);
figFRC = gcf;
%% Verification: Collocation using <https://sourceforge.net/p/cocotools/wiki/Home/ coco>
% Dankowicz, H., & Schilder, F. (2013).  _Recipes for Continuation,_ SIAM Philadelphia. 
% <https://doi.org/10.1137/1.9781611972573 <https://doi.org/10.1137/1.9781611972573>>

nCycles = 10;

coco = cocoWrapper(DS, nCycles, outdof);
set(coco,'initialGuess','forward')
set(coco.Options, 'NAdapt', 1);
set(coco.Options,'ItMX',10,'NTST', 30,'PtMX',200); %for convergence, smaller stepsize

figure(figFRC)
hold on;
startcoco = tic;
bd = coco.extract_FRC(OmegaRange);
timings.cocoFRC = toc(startcoco);