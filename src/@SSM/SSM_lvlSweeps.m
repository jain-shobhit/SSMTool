function FRC = SSM_lvlSweeps(obj, omRange, epSamp, ORDER)
% SSM_LVLSWEEPS This function performs a family of continuations of
% periodic orbits of slow dynamics. Continuation is performed with varied
% forcing frequency for each forcing amplitude.
% The continuation here starts from the guess of initial solution.
%
% FRC = SSM_LVLSWEEPS(OBJ,OMRANGE,EPSAMPS,ORDER)
%
% omRange:      continuation domain of forcing frequency, which should be near the
%               value of natural frequency with index 1
% epSamp :      sampled forcing amplitudes for forced response curves
% ORDER  :      expansion orders of SSM, can be scalar or array
%
% FRC:          FRC data struct
%
% See also: FRC_LEVEL_SET, SSM_POSWEEPS, SSM_EPSWEEPS



totalComputationTime = zeros(size(ORDER));

for j = 1:numel(ORDER)
    order = ORDER(j);
    startFRC = tic;
    if isempty(obj.System.spectrum)
        [~,~,~] = obj.System.linear_spectral_analysis();
    end
    lambda  = obj.System.spectrum.Lambda;
    assert(~isreal(lambda),'One or more eigenvalues must be underdamped for FRC computation using SSMs')
    
    % detect resonant eigenvalues in the parameter range
    if obj.System.order == 2
        assert(~isempty(obj.System.fext.epsilon), 'The epsilon field is empty in the dynamical system external forcing');
    else
        assert(~isempty(obj.System.Fext.epsilon), 'The epsilon field is empty in the dynamical system external forcing');
    end
    [resLambda,resFreq] = find_eigs_in_freq_range(omRange,lambda,obj.FRCOptions.resType);
    % obtain subintervals around each resonant eigenvalue
    [parNodes, nSubint] = subdivide_freq_range(omRange, resFreq);
    
    
    
    FRC = cell(nSubint,1);
    for i=1:nSubint
        % tune subinterval
        parSubRange = parNodes(i:i+1)';
        parSubRange = tune_parameter_range(parSubRange, obj.FRCOptions.frac, i, nSubint);
        
        % detect modes resonant with resLambda(i)
        [resModes,] = detect_resonant_modes(resLambda(i),lambda, obj.Options.IRtol);
        
        %% FRC computation within the subinterval
        disp('*****************************************');
        disp(['Calculating FRC using SSM with master subspace: [' num2str(resModes(:).') ']']);
        
        FRC{i} = lvlSweep(obj,resModes,order,parSubRange,epSamp);
        
        sweeps_plot(FRC{i},obj.FRCOptions.outdof,order);
    end
end
% concatenate cell contents as struct arrays
FRC = cat(1,FRC{:});

totalComputationTime(j) = toc(startFRC);

for j = 1:numel(ORDER)
    disp(['Total time spent on FRC computation upto O(' num2str(ORDER(j)) ') = ' datestr(datenum(0,0,0,0,0,totalComputationTime(j)),'HH:MM:SS')])
end
end

function parRange = tune_parameter_range(parRange, frac, i, nSubint)
% amplify the parameter subinterval except on the first and the last nodes.
if i>1
    parRange(1) = frac(1)*parRange(1);
end
if i<nSubint
    parRange(2) = frac(2)*parRange(2);
end
end

function [resLambda, resFreq] = find_eigs_in_freq_range(Omega,lambda,resType)
switch resType
    case '1:1'
        natFreq = imag(lambda);
        resFreqID = intersect(find(natFreq>Omega(1)), find(natFreq<Omega(end)));
        resFreq = natFreq(resFreqID);
        resLambda = lambda(resFreqID);
    case '2:1'  % Case of subharmonic resonance
        
        natFreq = imag(lambda);
        resFreqID = intersect(find(2*natFreq>Omega(1)), find(2*natFreq<Omega(end)));
        resFreq = natFreq(resFreqID);
        resLambda = lambda(resFreqID);
end
% remove repetitive eigenvalues, e.g., 1:1 internal resonance
dFreq = resFreq(2:end)-resFreq(1:end-1);
idrep = abs(dFreq)<1e-3*resFreq(1:end-1);
resFreq(idrep)   = [];
resLambda(idrep) = [];
assert(~isempty(resFreq),'Input frequency range should include at least one  (multiple) natural frequency'); % we could still program this case
end

function [freqNodes, nSubint] = subdivide_freq_range(parRange,resFreq)
nSubint = numel(resFreq);
freqNodes = 0.5*(resFreq(1:end-1)+resFreq(2:end)); % center points of two adjacent resonant modes
freqNodes = [parRange(1); freqNodes; parRange(2)];
end

function [FRC] = lvlSweep(obj, resModes, order, parRange, epsSamp)
%  EXTRACT_FRC This function extracts the steady-state amplitude, phase shift and frequency
%
% for periodically forced systems from *Forced response curves*. Current implementation
% assumes forcing of the form $\mathbf{F}_0\cos(m\Omega t) = \mathbf{F}_m^{ext}e^{\iota
% m\Omega t}+\bar{\mathbf{F}}_m^{ext}e^{-\iota m\Omega t}$. For two-dimensional
% SSMs, we use the normal form of paramaterization, where we choose the following
% form of autonomous reduced dynamics as
%
% $$\mathbf{R}_{0}(\mathbf{p})=\left[\begin{array}{c}\lambda p\\\bar{\lambda}\bar{p}\end{array}\right]+\sum_{j=1}^{M}\left[\begin{array}{c}\gamma_{j}p^{j+1}\bar{p}^{j}\\\bar{\gamma}_{j}p^{j}\bar{p}^{j+1}\end{array}\right],$$
%
% and assuming that the $m$-th harmonic of forcing frequency is nearly resonant
% with the modal subspace eigenvalues, i.e., $\lambda-\iota m\Omega,\bar{\lambda}+\iota
% m\Omega\approx0$. Denoting $\mathbf{W}_{\mathcal{M}}=[\mathbf{w},\,\bar{\mathbf{w}}]\in\mathbb{C}^{N\times2}$,
% we obtain the leading-order, non-autonomous reduced dynamics using the normal
% form parametrization as
%
% $$\mathbf{R}_{1,0}(\phi)=\left[\begin{array}{c}\mathbf{w}^{\star}\mathbf{F}_{m}^{ext}e^{\iota
% m\phi}\\\mathbf{\bar{w}}^{\star}\bar{\mathbf{F}}_{m}^{ext}e^{-\iotam\phi}\end{array}\right],\quad
% \dot{\phi}=\Omega$$
%
% We are able to explicitly express the *forced response curve* as the zero
% level set of the functions
%
% $$\mathcal{G}^{\pm}(\rho,\Omega,|f|):=b(\rho,m\Omega)\rho\pm\sqrt{|f|^{2}-\left[a(\rho)\right]^{2}},$$
%
% where
%
% $$a(\rho)=\sum_{j=1}^{M}\Re(\gamma_{j})\rho^{2j+1}+\rho\Re(\lambda),$$
%
% $$b(\rho,m\Omega)=\sum_{j=1}^{M}\Im(\gamma_{j})\rho^{2j}+\Im(\lambda)-m\Omega$,
%
% $$$f	=\epsilon\frac{\mathbf{w}^{\star}\mathbf{F}_{0}}{2}.$$
%
dimModes = numel(resModes);
assert(dimModes==2,'levelset method for FRC extraction is valid only for 2-dimensional SSMs. Please use the continuation method')

% compute two dimensional autonomous SSM
obj.choose_E(resModes)
[W0, R0] = obj.compute_whisker(order);

% autonomous reduced dynamics coefficients
gamma  = compute_gamma(R0);
lambda = obj.E.spectrum(1);

% get options
[nt, nRho, nPar, rhoScale, nPsi, outdof, saveIC]  = ...
    deal(obj.FRCOptions.nt, obj.FRCOptions.nRho, ...
    obj.FRCOptions.nPar, obj.FRCOptions.rhoScale, obj.FRCOptions.nPsi,...
    obj.FRCOptions.outdof, obj.FRCOptions.saveIC);

% Initialize FRC as a struct array
par = linspace(parRange(1),parRange(end), nPar);
nSamp = numel(epsSamp);
FRC = cell(nPar,nSamp);


%% Grid for evaluation of $\mathbf{r}(rho,psi,Omega)$


Omega_init   = par(1);
omegaRange   = par;

% leading-order modal forcing coefficient
f                    = compute_f(obj,W0,R0,0,Omega_init);

[rho, psi, RHO, PSI] = compute_polar_grid(omegaRange,epsSamp,gamma,lambda,f,rhoScale, nRho, nPsi);

%% Obtain FRC

parfor j = 1:nPar
    
    FRCj = FRC{j,:};
    
    Omega   = par(j);
    
    
    % compute non-autonomous SSM coefficients
    [W1, R1] = obj.compute_perturbed_whisker(order-1,W0,R0,Omega);
    
    ii = 1;
    for epsilon = epsSamp
        % reduced dynamics on the grid
        [rhodot, rhopsidot, eta] = compute_reduced_dynamics_2D_polar(RHO{ii},PSI{ii}, ...
            lambda, gamma, R1,Omega,epsilon);
        
        % Numerical computation of fixed points of the reduced dynamics
        [rho0, psi0] = compute_fixed_points_2D(rho{ii}, psi{ii}, rhodot, rhopsidot);
        
        % Stabilty calculation
        stability    = check_stability(rho0,psi0,gamma,lambda,epsilon,R1);
        
        % output data structure
        FRCj{ii} = compute_output_polar2D(rho0,psi0,stability,epsilon,Omega*ones(size(rho0)),W0,W1,eta,nt, saveIC, outdof);
        ii = ii+1;
    end
    
    FRC(j,:) = FRCj;
end

FRC = cat(1,FRC{:});
end

function [rho, psi, RHO, PSI] = compute_polar_grid(omegaRange,epsSamp,gamma,lambda,f,rhoScale, nRho, nPsi)
ii = 1;
for epsilon = epsSamp
    % *Explicit quadratic approximation of the backbone curve*
    %
    % $$\rho_{backbone} = \sqrt{\frac{\Omega-\Im(\lambda)}{\Im(\gamma_1)}}$$
    rhomaxBB = max(real(sqrt((omegaRange - imag(lambda))/imag(gamma(1)))));
    rhomaxlin = full(abs(epsilon * f/real(lambda)));
    
    if isempty(f)
        assert( rhomaxBB ~= 0, 'Estimation of maximal amplitude failed')
        rhomax = rhoScale * rhomaxBB;
        
    else
        if rhomaxBB ~=0
            % heuristically choose maximum value of polar radius as follows
            rhomax = rhoScale * min(rhomaxBB, rhomaxlin);
        else
            rhomax = rhoScale * rhomaxlin;
        end
    end
    %compute grids
    rhoi = (rhomax/nRho) * (1:nRho);
    psii = (2*pi/nPsi) * (0:nPsi-1);
    [RHOi,PSIi] = meshgrid(rhoi,psii);
    
    rho{ii} = rhoi;
    psi{ii} = psii;
    RHO{ii} = RHOi;
    PSI{ii} = PSIi;
    ii = ii+1;
end
end

function [f] = compute_f(obj,W0,R0,order,Omega)
% compute non-autonomous SSM coefficients
[~, R1] = obj.compute_perturbed_whisker(order-1,W0,R0,Omega);
f =  nonzeros(R1(1).R(1).coeffs); % leading-order modal forcing coefficient
end

function sweeps_plot(FRCom,outdof,order)

omega = [FRCom.Omega];
epsilon = [FRCom.epsilon];
Aout = [FRCom.Aout];
stab = [FRCom.stability];
Znorm = [FRCom.Znorm];

numOutdof = numel(outdof);
numPts    = numel(stab);
Aout = reshape(Aout, [numOutdof, numPts]);
Aout = Aout';

plot3_frc_full(omega,epsilon,Znorm,outdof,Aout,stab,order)

end
