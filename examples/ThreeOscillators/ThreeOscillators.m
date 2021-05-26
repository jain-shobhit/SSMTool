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
% *1:1:1 internal resonance*

S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
order = 3;
outdof = [1 2 3];
%% 
% We are interested in the FRC over the frequency span [0.995 1.016]. We can 
% call the routine extract_FRC to extract the FRC. All the three natural frequencies 
% are inside the frequency span. The frequency span is divided into a single subinterval 
% because of the 1:1:1 internal resonance. A single continuation run is involved 
% to get the FRC for the subinterval.

set(S.Options, 'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'method','continuation ep')
set(S.FRCOptions, 'outdof',outdof)
set(S.FRCOptions, 'initialSolver', 'fsolve');
set(S.contOptions,'PtMX',250);

freqRange = [0.995 1.016];
% call extract_FRC to calculate the FRC
FRC = S.extract_FRC('freq',freqRange,order);
%% 
% On the other hand, we can use SSM-ep toolbox to obtaint the FRC with a single 
% continuation run as well given $\Omega\approx1$ in the whole frequency span. 
% The SSM-ep toolbox is the minic of the ep-toolbox of COCO. It is in the middle 
% level of the SSMTool while the extract_FRC routine is in the high level of the 
% SSMTool. The continuation-based method implemented in the extract_FRC routine 
% is actually built on the SSM-ep toolbox. Users may use the extract_FRC routine 
% initially. They may explore the SSM-ep toolbox given this toolbox has more functionalities, 
% e.g., bifurcation analysis.
% 
% We call the function SSM_isol2ep, whose arguments can be found via help

help SSM_isol2ep
%% 
% In this example, the resonant_modes should be [1 2 3 4 5 6] due to all modes 
% are involved in the resonance. We have $\omega_1=\omega_2=\omega_3\approx\Omega$ 
% and then $\mathbf{r}=[1,1,1]$ (mfreqs).  The argument parName can be amp or 
% freq. When parName='freq'/'amp', the forced response curve with varied excitation 
% frequency $\Omega$ / excitation amplitude $\epsilon$ is obtained.

resonant_modes = [1 2 3 4 5 6];
mFreq  = [1 1 1];
set(S.contOptions,'PtMX',250);
startep = tic;
FRC_ep_polar = S.SSM_isol2ep('isol_polar',resonant_modes, order, mFreq, 'freq', freqRange, outdof);
timings.epPolarFRC = toc(startep);
%% 
% As seen in the continuation hisotry, the continuation run terminated at a 
% point where $\rho_3\to0$ (denoted as MX point), which triggers the singularity 
% of the vector field. As an alternative, we can use Cartesian coordinates to 
% remove the MX point.

set(S.FRCOptions, 'coordinates','cartesian');
startep = tic;
FRC_ep_cart = S.SSM_isol2ep('isol_cart',resonant_modes, order, mFreq, 'freq', freqRange, outdof);
timings.epCartFRC = toc(startep);
%% Validation using COCO
% To conclude this example, we use the po-toolbox (collocation method) of COCO 
% to validate the results obtained from the SSM analysis. As we can see, the results 
% of the two methods match well. In addition, the runtime of collocation method 
% is nearly three times of that of the SSM analysis.

nCycles = 100;
coco = cocoWrapper(DS, nCycles, outdof);
set(coco.Options, 'PtMX', 500);
set(coco.Options, 'NAdapt', 1, 'h_max', 50);

startcoco = tic;
bd2 = coco.extract_FRC(freqRange);
timings.cocoFRC = toc(startcoco)