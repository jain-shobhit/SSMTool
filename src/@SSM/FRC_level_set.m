function [FRC] = FRC_level_set(obj, resModes, order, parName, parRange)
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
gamma = compute_gamma(R0);
lambda = obj.E.spectrum(1);

% get options
[nt, nRho, nPar, rhoScale, nPsi, outdof, saveIC]  = ...
    deal(obj.FRCOptions.nt, obj.FRCOptions.nRho, ...
    obj.FRCOptions.nPar, obj.FRCOptions.rhoScale, obj.FRCOptions.nPsi,...
    obj.FRCOptions.outdof, obj.FRCOptions.saveIC);

% Initialize FRC as a struct array
par = linspace(parRange(1),parRange(end), nPar);

FRC = cell(nPar,1);
for j = 1:nPar
    % setup
    switch lower(parName)
        case 'freq' % vary Omega keeping epsilon constant
            Omega = par(j);
            epsilon = obj.System.fext.epsilon;
            obj.System.Omega = Omega;
        case 'amp' % vary epsilon keeping Omega constant
            Omega = obj.System.Omega;
            epsilon = par(j);
            obj.System.fext.epsilon = epsilon;
    end
    if obj.System.Options.BaseExcitation
        epsilon = epsilon*Omega^2;
    end
    % compute non-autonomous SSM coefficients
    [W1, R1] = obj.compute_perturbed_whisker(order);
    
    %     compute grid for rho and psi
    f =  R1{1}.coeffs(1,2); % leading-order modal forcing coefficient    
    
    [rho, psi, RHO, PSI] = compute_polar_grid(Omega,epsilon,gamma,lambda,f,rhoScale, nRho, nPsi);
    
    % reduced dynamics on the grid
    [rhodot, rhopsidot, eta] = compute_reduced_dynamics_2D_polar(RHO,PSI, ...
        lambda, gamma, R1,Omega,epsilon);
    
    % Numerical computation of fixed points of the reduced dynamics
    [rho0, psi0] = compute_fixed_points_2D(rho, psi, rhodot, rhopsidot);
    
    % Stabilty calculation
    stability = check_stability(rho0,psi0,gamma,lambda,epsilon,f);
    
    
    % output data structure
    FRC{j} = compute_output_polar2D(rho0,psi0,stability,epsilon,Omega*ones(size(rho0)),W0,W1,eta,nt, saveIC, outdof);    
end

FRC = cat(1,FRC{:});
end
%%
function [rho, psi, RHO, PSI] = compute_polar_grid(Omega,epsilon,gamma,lambda,f,rhoScale, nRho, nPsi)
%%
% *Explicit quadratic approximation of the backbone curve*
%
% $$\rho_{backbone} = \sqrt{\frac{\Omega-\Im(\lambda)}{\Im(\gamma_1)}}$$
rhomaxBB = real(sqrt((Omega - imag(lambda))/imag(gamma(1))));
rhomaxlin = full(abs(epsilon * f/real(lambda)));
%
if rhomaxBB ~=0
    % heuristically choose maximum value of polar radius as follows
    rhomax = rhoScale * min(rhomaxBB, rhomaxlin);
else
    rhomax = rhoScale * rhomaxlin;
end
%compute grids
rho = (rhomax/nRho) * (1:nRho);
psi = (2*pi/nPsi) * (0:nPsi-1);
[RHO,PSI] = meshgrid(rho,psi);
end

