function varargout = extract_ridges_trenches(obj,oid,resonant_modes,order,mFreqs,parRange,outdof,optdof,varargin)
% DAMPED_BACKBONE_CURVE_CONT_EP This function extract the damped backbone
% curve (defined as amplitude resonance) of reduced-order model (in slow
% phase) as optimization problem of fixed points. Note that this backbone
% curve will be slightly deviated from the one for a given physical
% coordinate because the contribution of non-autnomous SSM to the response
% in physical coordinates. Each equilibirum point corresponds to a periodic
% orbit in full system. The continuation here starts from the guess of
% initial solution (for fixed point, adjoints are initialization-free).
%
% FRC = DAMPED_BACKBONE_CURVE_CONT_EP(OBJ,OID,RESONANT_MODES,ORDER,MFREQS,PARRANGE,OUTDOF,OPTDOF,VARARGIN)
%
% oid:      runid of continuation
% resonant_modes:    master subspace
% order:    expansion order of SSM
% mFreqs:   internal resonance relation vector
% optobj:   [] or a structure with fields {optfunc, doptfunc, name}, where
%           optfunc is the optimization objective function
%           doptfunc is the Jacobian of objective function, and
%           name: a string for the fuction name that used in plotting
% parRange: continuation domain of parameters {[freq1,freq2],[eps1,eps2]},
%           the frequency range should be near the
%           value of natural frequency with index 1 in the mFreq
% outdof:   output for dof in physical domain
% varargin: [{p0,z0}], ['saveICs'] where {p0,z0} are initial solution
%           guesses and saveICs is a flag saving a point on trajectory as initial
%           condition for numerical integration


m = numel(mFreqs);
assert(numel(resonant_modes)==2*m, 'The master subspace is not %dD.',2*m);
nCycle = obj.FRCOptions.nCycle;
isl2norm = strcmp(obj.FRCOptions.DBCobjnorm, 'l2');
if ~isl2norm
   assert(numel(optdof)==1,'Number of optdof should be one when Linfty norm is used'); 
end
%% Checking whether internal resonance indeed happens
if isempty(obj.System.spectrum)
    [~,~,~] = obj.System.linear_spectral_analysis();
end
% eigenvalues Lambda is sorted in descending order of real parts
% positive imaginary part is placed first in each complex pair
lambda = obj.System.spectrum.Lambda(resonant_modes);
lambdaRe = real(lambda);
lambdaIm = imag(lambda);
check_spectrum_and_internal_resonance(lambdaRe,lambdaIm,mFreqs);

%% SSM computation of autonomous part
obj.choose_E(resonant_modes)
% compute autonomous SSM coefficients
[W_0,R_0] = obj.compute_whisker(order);

% check reduced dynamics (consistency with internal resonance)
[beta,kappa] = check_auto_reduced_dynamics(R_0,order,mFreqs);

% assemble expansion coefficients and indices
wcoeffs = []; cind = []; dind = [];
for k=1:order
    W = W_0(k);
    coeffs = W.coeffs(optdof,:);
    ind = W.ind;
    if ~isempty(coeffs)
        wcoeffs = [wcoeffs coeffs];
        cind = [cind; ind(:,1:2:end-1)];
        dind = [dind; ind(:,2:2:end)];
    end
end

%% SSM computation of non-autonomous part
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
kNonauto = []; % (pos) kappa indices with resonance
kappa_set= [obj.System.Fext.data.kappa];
kappa_pos = kappa_set(kappa_set>0);
num_kappa = numel(kappa_pos); % number of kappa pairs
for k=1:num_kappa
    kappak = kappa_pos(k);
    idm = find(mFreqs(:)==kappak); % idm could be vector if there are two frequencies are the same
    obj.System.Omega = lambdaIm(2*idm(1)-1);
    
    [W_1, R_1, Flead] = obj.compute_perturbed_whisker(0,[],[]);

    idk = find(kappa_set==kappak);
    R_10 = R_1(idk).R.coeffs;
    r = R_10(2*idm-1);
    
    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];  
    kNonauto = [kNonauto; idk];
end

%% Construct COCO-compatible vector field 
lamd  = struct();
lamd.lambdaRe = lambdaRe; lamd.lambdaIm = lambdaIm;
Nonauto = struct();
Nonauto.iNonauto = iNonauto; Nonauto.rNonauto = rNonauto; Nonauto.kNonauto = kNonauto;
[fdata,data_dir] = create_reduced_dynamics_data(beta,kappa,lamd,mFreqs,Nonauto,W_0,W_1,order,resonant_modes);
ispolar = strcmp(obj.FRCOptions.coordinates, 'polar');
fdata.ispolar = ispolar;
fdata.isbaseForce = obj.System.Options.BaseExcitation;
if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
    odefun_dfdx = @(z,p) ode_2mDSSM_polar_DFDX(z,p,fdata);
    odefun_dfdp = @(z,p) ode_2mDSSM_polar_DFDP(z,p,fdata);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
    odefun_dfdx = @(z,p) ode_2mDSSM_cartesian_DFDX(z,p,fdata);
    odefun_dfdp = @(z,p) ode_2mDSSM_cartesian_DFDP(z,p,fdata);    
end
funcs  = {odefun,odefun_dfdx,odefun_dfdp};

%% create data structure for optimization objective
optData = create_data_for_po_amp(cind,dind,mFreqs,wcoeffs,optdof,obj,Flead,isl2norm);

%% continuation of reduced dynamics w.r.t. parName
% construct initial solution
if obj.System.order==2
    p0 = [obj.System.Omega; obj.System.fext.epsilon];
else
    p0 = [obj.System.Omega; obj.System.Fext.epsilon];
end
[p0,z0] = initial_fixed_point(p0,obj.FRCOptions.initialSolver,ispolar,...
    odefun,nCycle,m,varargin{:});
if strcmp(obj.FRCOptions.DBCobjnorm, 'l2')
    [runid,runidbc,args1,args2] = damped_backbone_l2(oid,obj.contOptions,...
        obj.FRCOptions,funcs,z0,p0,ispolar,optData,parRange);
else
    [runid,runidbc,args1,args2] = damped_backbone_linf(oid,obj.contOptions,...
        obj.FRCOptions,funcs,z0,p0,ispolar,optData,parRange);    
end
%% plot dbc (in reduced coordinates)
numBP = plot_dbc_figures(runidbc,runid);

%% extract results of reduced dynamics at sampled frequencies
FRCdata = struct();        FRCdata.isomega = true;
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order;
FRCdata.ispolar = ispolar; FRCdata.modes  = resonant_modes;
FRCs = cell(numBP,1);
for idbp = 1:numBP
    FRC = ep_reduced_results(runidbc{idbp},obj.FRCOptions.sampStyle,ispolar,true,args1,args2,'no-plot');
    FRC = FRC_reduced_to_full(obj,Nonauto,'ep',FRC,FRCdata,W_0,W_1,outdof,varargin{:});
    FRCs{idbp} = FRC;
end

varargout{1} = FRCs;
fdir = fullfile(data_dir,runid,'SSMep.mat');
save(fdir, 'FRCs');
end

function numBP = plot_dbc_figures(runidbc,runid)

thm = struct( 'special', {{'BP'}});
thm.lspec = {{'r-'  'LineWidth'  [1]},{'r--'  'LineWidth'  [1]}};

figure;
coco_plot_bd(thm,runid,'om','eps','obj', @(x) abs(x)); hold on
numBP = numel(runidbc);
for idbp = 1:numBP
    coco_plot_bd(runidbc{idbp},'om','eps','obj', @(x) abs(x)); hold on
end
grid on
xlabel('$\Omega$','interpreter','latex','FontSize',16);
ylabel('$\epsilon$','interpreter','latex','FontSize',16);
zlabel('$||z_\mathrm{opt}||$','interpreter','latex','FontSize',16);

end