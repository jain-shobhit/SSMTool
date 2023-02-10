function BB = extract_backbone(obj, modes, omegaRange, order, varargin)
%  EXTRACT_BACKBONE This function extracts the *Backbone curves in Polar coordinates.* For two-dimensional
% SSMs, we use the normal form of paramaterization, where we choose the following
% form of autonomous reduced dynamics as
%
% $$\mathbf{R}_{0}(\mathbf{p})=\left[\begin{array}{c}\lambda p\\\bar{\lambda}\bar{p}\end{array}\right]+\sum_{j=1}^{M}\left[\begin{array}{c}\gamma_{j}p^{j+1}\bar{p}^{j}\\\bar{\gamma}_{j}p^{j}\bar{p}^{j+1}\end{array}\right],$$
%
% Subsitution of $p=\rho e^{\mathrm{i}\theta}$ and $\bar{p}=\rho e^{-\mathrm{i}\theta}$
% to $\dot{\mathbf{p}}=\mathbf{R}_0(\mathbf{p})$ yields
%
% $$\dot{\rho}=a(\rho), \dot{\theta}=b(\rho)$$
%
% where
%
% $$a(\rho)=\sum_{j=1}^{M}\Re(\gamma_{j})\rho^{2j+1}+\rho\Re(\lambda),$$
%
% $$b(\rho)=\sum_{j=1}^{M}\Im(\gamma_{j})\rho^{2j+1}+\rho\Im(\lambda).$$
%
% It follows that the _backbone curves_ in polar coordinates is given by $\Omega=\frac{b(\rho)}{\rho}$.
%
% The range of rho is determined by quadratic approximation of backbone
% curve if varargin is empty. Otherwise, it is specified via varargin

f1 = figure('Name','Norm');
if isnumeric(obj.FRCOptions.outdof)
    f2 = figure('Name',['Amplitude at DOFs ' num2str(obj.FRCOptions.outdof(:)')]);
else
    f2 = figure('Name','Amplitude at DOFs');
end
figs = [f1, f2];
colors = get(0,'defaultaxescolororder');


assert(numel(modes)==2,'The analytic backbone computation can only be performed for a two-dimensional SSM/LSM')
% get options
[nt, nRho, nOmega, rhoScale, outdof, saveIC]  = ...
    deal(obj.FRCOptions.nt, obj.FRCOptions.nRho, ...
    obj.FRCOptions.nPar, obj.FRCOptions.rhoScale, ...
    obj.FRCOptions.outdof, obj.FRCOptions.saveIC);

%% setup
startBB = tic;
obj.choose_E(modes);
lambda = obj.E.spectrum(1);

% some checks
assert(~isreal(lambda),'The eigenvalues associated to the modal subspace must be complex for analytic backbone computation')
omega0 = abs(imag(lambda));
assert(prod([omega0-omegaRange(1),omega0-omegaRange(end)])<0,'The supplied omegaRange must contain the natural frequency associated to the modes')

%% loop over orders
norders = numel(order);
for k=1:norders
    %% compute autonomous SSM coefficients
    [W0,R0] = obj.compute_whisker(order(k));
    gamma = compute_gamma(R0);

    %% compute backbone
    if numel(varargin)==0
        rho = compute_rho_grid(omegaRange,nOmega,rhoScale,gamma,lambda,nRho);
    else
        rhomax = varargin{1};
        rho = linspace(0.001*rhomax,rhomax,nRho);
    end
    [~,b] = frc_ab(rho, 0, gamma, lambda);
    omega = b./rho;
    idx = [find(omega<omegaRange(1)) find(omega>omegaRange(2))];
    rho(idx) = []; omega(idx) = [];

    %% Backbone curves in Physical Coordinates
    stability = true(size(rho)); psi = zeros(size(rho)); epsilon = 0;
    BB = compute_output_polar2D(rho,psi,stability,epsilon,omega,W0,[],1,nt, saveIC, outdof);

    %% plotting
    plot_FRC(BB,outdof,order(k),'freq','lines',figs,colors(k,:));
end
totalComputationTime = toc(startBB);
disp(['Total time spent on backbone curve computation = ' datestr(datenum(0,0,0,0,0,totalComputationTime),'HH:MM:SS')])
end

function rho = compute_rho_grid(omegaRange,nOmega,rhoScale,gamma,lambda,nRho)
omega = linspace(omegaRange(1),omegaRange(end), nOmega);
%%
% *Explicit quadratic approximation of the backbone curve*
%
% $$\rho_{backbone} = \sqrt{\frac{\Omega-\Im(\lambda)}{\Im(\gamma_1)}}$$
rho_bb = real(sqrt((omega - imag(lambda))/imag(gamma(1))));
rhomax = rhoScale * max(rho_bb);
rho = (rhomax/nRho) * (1:nRho);
end