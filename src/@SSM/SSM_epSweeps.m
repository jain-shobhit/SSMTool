function varargout = SSM_epSweeps(obj,oid,resonant_modes,order,mFreqs,epSamp,omRange,outdof,varargin)
% SSM_EPSWEEPS This function performs a family of continuation of 
% equilibrium points of slow dynamics. In the first continuation run,
% forcing amplitude is changed and equilibria at sampled forcing
% frequencies will be lablled. Then continuation is performed with varied
% forcing frequency starting from each sampled equilibria. Note that
% equilibrium points in slow dynamics corresponds to a periodic orbit in
% the regular time dynamics. The continuation here starts from the guess of
% initial solution. 
%
% FRC = SSM_EPSWEEPS(OBJ,OID,MODES,ORDER,MFREQS,EPSAMPS,OMRANGE,OUTDOF,VARARGIN)
%
% oid:      runid of continuation
% modes:    master spectral subspace
% order:    expansion order of SSM
% mFreqs:   internal resonance relation vector
% epSamp:   sampled forcing amplitudes for forced response curves
% omRange:  continuation domain of forcing frequency, which should be near the
%           value of natural frequency with index 1
% outdof:   output for dof in physical domain
% varargin: [{p0,z0}], ['saveICs'] where {p0,z0} are initial solution
%           guesses and saveICs is a flag saving a point on trajectory as initial
%           condition for numerical integration

m = numel(mFreqs);
assert(numel(resonant_modes)==2*m, 'The master subspace is not %dD.',2*m);
nOmega = obj.FRCOptions.nPar;
nCycle = obj.FRCOptions.nCycle;
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

% check reduced dynamics to see its consistent with reduced dynamics
[beta,kappa] = check_auto_reduced_dynamics(R_0,order,mFreqs);

%% SSM computation of non-autonomous part
iNonauto = []; % indices for resonant happens
rNonauto = []; % value of leading order contribution
kNonauto = []; % (pos) kappa indices with resonance
kappa_set= obj.System.Fext.kappas; % each row corresponds to one kappa
kappa_pos = kappa_set(kappa_set>0);
num_kappa = numel(kappa_pos); % number of kappa pairs
for k=1:num_kappa
    kappak = kappa_pos(k);
    idm = find(mFreqs(:)==kappak); % idm could be vector if there are two frequencies are the same
    obj.System.Omega = lambdaIm(2*idm(1)-1);
    
    [W_1, R_1] = obj.compute_perturbed_whisker(order);

    R_10 = R_1{1}.coeffs;
    idk = find(kappa_set==kappak);
    r = R_10(2*idm-1,idk);
    
    iNonauto = [iNonauto; idm];
    rNonauto = [rNonauto; r];  
    kNonauto = [kNonauto; idk];
end
%% Construct COCO-compatible vector field 
lamd  = struct(); 
lamd.lambdaRe = lambdaRe; lamd.lambdaIm = lambdaIm;
Nonauto = struct();
Nonauto.iNonauto = iNonauto; Nonauto.rNonauto = rNonauto; Nonauto.kNonauto = kNonauto;
[fdata,~] = create_reduced_dynamics_data(beta,kappa,lamd,mFreqs,Nonauto,W_0,W_1,order,resonant_modes);

ispolar = strcmp(obj.FRCOptions.coordinates, 'polar');
fdata.ispolar = ispolar;
fdata.isbaseForce = obj.System.Options.BaseExcitation;
if ispolar
    odefun = @(z,p) ode_2mDSSM_polar(z,p,fdata);
else
    odefun = @(z,p) ode_2mDSSM_cartesian(z,p,fdata);
end
funcs  = {odefun};
% 
%% continuation of reduced dynamics w.r.t. epsilon
prob = coco_prob();
prob = cocoSet(obj.contOptions, prob);
% construct initial solution
if obj.System.order==2
    p0 = [obj.System.Omega; obj.System.fext.epsilon];
else
    p0 = [obj.System.Omega; obj.System.Fext.epsilon];
end
[p0,z0] = initial_fixed_point(p0,obj.FRCOptions.initialSolver,ispolar,...
    odefun,nCycle,m,varargin{:});
% call ep-toolbox
prob = ode_isol2ep(prob, '', funcs{:}, z0, {'om' 'eps'}, p0);
% define monitor functions to state variables
[prob, args1, args2] = monitor_states(prob, ispolar, m);
prob  = coco_add_event(prob, 'UZ', 'eps', epSamp); % label epSamp with UZ
runid = [oid, 'eps'];
runid = coco_get_id(runid, 'ep');
fprintf('\n Run=''%s'': Continue equilibria with varied epsilon.\n', ...
  runid);
cont_args = [{'eps'},args1(:)' ,args2(:)',{'om'}];
coco(prob, runid, [], 1, cont_args, [0.9*min(epSamp) 1.1*max(epSamp)]);

%% a family of continuation of reduced dynamics w.r.t. Omega
bd = coco_bd_read(runid);
sampLabs = coco_bd_labs(bd, 'UZ'); 
numLabs  = numel(sampLabs);
if numLabs~=numel(epSamp)
    warning('The number of sampled forcing amplitude is not the same as requested.');
end
ep   = coco_bd_vals(bd, sampLabs, 'eps');
Labs  = zeros(numLabs,1);

for k=1:numLabs
    epk     = epSamp(k);
    [~,id]  = min(abs(ep-epk));
    Labs(k) = id;
    prob = coco_prob();
    prob = cocoSet(obj.contOptions,prob);
    prob = ode_ep2ep(prob, '', runid, sampLabs(id));
    
    if strcmp(obj.FRCOptions.sampStyle, 'uniform')
        omSamp = linspace(omRange(1),omRange(2), nOmega);
        prob   = coco_add_event(prob, 'UZ', 'om', omSamp);
    end
    
    [prob, args1, args2] = monitor_states(prob, ispolar, m);
    cont_args = [{'om'},args1(:)' ,args2(:)',{'eps'}];
    runidk = [oid,'eps',num2str(k)];
    runidk = coco_get_id(runidk, 'ep');

    fprintf('\n Run=''%s'': Continue equilibria with varied omega at eps equal to %d.\n', ...
      runidk, ep(id));
  
    coco(prob, runidk, [], 1, cont_args, omRange);
end        

%% results extraction
FRCom = cell(numLabs,1);
FRCdata = struct();        FRCdata.isomega = true; 
FRCdata.mFreqs  = mFreqs;  FRCdata.order  = order; 
FRCdata.ispolar = ispolar; FRCdata.modes  = resonant_modes;
for k=1:numLabs
    %% results in reduced system
    runidk = [oid,'eps',num2str(k)];
    runidk = coco_get_id(runidk, 'ep');
    FRC = ep_reduced_results(runidk,obj.FRCOptions.sampStyle,ispolar,true,args1,args2,'plot-off');
    %% results in physical domain
    if ~isempty(outdof)
    fprintf('Calculate FRC in physical domain at epsilon %d\n', ep(Labs(k)));
    FRC = FRC_reduced_to_full(obj,Nonauto,'ep',FRC,FRCdata,W_0,W_1,outdof,varargin{:});
    end
    FRCom{k} = FRC;
end

% save results
FRC = struct();
FRC.SSMorder   = order;
FRC.SSMnonAuto = obj.Options.contribNonAuto;
FRC.SSMispolar = ispolar;
FRC.FRCom = FRCom;
FRC.eps = ep(Labs);
varargout{1} = FRC;
fdir = fullfile(pwd,'data',runid,'SSMep.mat');
save(fdir, 'FRC');

%% Visualization of results
sweeps_plot(m,ispolar,numLabs,oid,FRCom,outdof,order);
end




function sweeps_plot(m,ispolar,numLabs,oid,FRCom,outdof,order)
% Plot FRC in normal coordinates
thm = struct( 'special', {{'SN', 'HB'}});
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};

for i=1:m
    % 2D plot
    figure;
    if ispolar
        rhoi = strcat('rho',num2str(i));
        for k=1:numLabs
            runidk = [oid,'eps',num2str(k)];
            runidk = coco_get_id(runidk, 'ep');
            coco_plot_bd(thm, runidk, 'om', rhoi, @(x) abs(x)); hold on
        end
    else
        Rezi = strcat('Rez',num2str(i));
        Imzi = strcat('Imz',num2str(i));
        for k=1:numLabs
            runidk = [oid,'eps',num2str(k)];
            runidk = coco_get_id(runidk, 'ep');
            coco_plot_bd(thm, runidk, 'om', {Rezi,Imzi}, @(x,y) sqrt(x.^2+y.^2)); hold on
        end              
    end
    grid on; box on; axis tight
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel(strcat('$\rho_',num2str(i),'$'),'interpreter','latex','FontSize',16);
    % 3D plot
    figure;
    if ispolar
        rhoi = strcat('rho',num2str(i));
        for k=1:numLabs
            runidk = [oid,'eps',num2str(k)];
            runidk = coco_get_id(runidk, 'ep');
            coco_plot_bd(thm, runidk, 'om', 'eps', rhoi, @(x) abs(x)); hold on
        end
    else
        Rezi = strcat('Rez',num2str(i));
        Imzi = strcat('Imz',num2str(i));
        for k=1:numLabs
            runidk = [oid,'eps',num2str(k)];
            runidk = coco_get_id(runidk, 'ep');
            coco_plot_bd(thm, runidk, 'om', 'eps', {Rezi,Imzi}, @(x,y) sqrt(x.^2+y.^2)); hold on
        end              
    end
    grid on; box on; axis tight
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\epsilon$','interpreter','latex','FontSize',16);
    zlabel(strcat('$\rho_',num2str(i),'$'),'interpreter','latex','FontSize',16);    
end

% Plot FRC in system coordinates
% setup legend
legTypes = zeros(numLabs,1); % each entry can be 1/2/3 - all stab/all unstab/mix
for k=1:numLabs
   stab = FRCom{k}.st;
   if all(stab)
       legTypes(k) = 1;
   elseif all(~stab)
       legTypes(k) = 2;
   else
       legTypes(k) = 3;
   end
end
legDisp = []; % figure with legends displayed
idmix = find(legTypes==3);
if ~isempty(idmix)
    legDisp = [legDisp idmix(1)];
else
    idallstab = find(legTypes==1);
    if ~isempty(idallstab)
        legDisp = [legDisp,idallstab(1)];
    end
    idallunstab = find(legTypes==2);
    if ~isempty(idallunstab)
        legDisp = [legDisp,idallunstab(1)];
    end
end
if ~isempty(outdof)
ndof = numel(outdof);
for i=1:ndof
    % 2D plot
    figure; hold on
    for k=1:numLabs
        om = FRCom{k}.om;
        Aout = FRCom{k}.Aout_frc(:,i);
        stab = FRCom{k}.st;
        if any(legDisp==k)
            stab_plot(om,Aout,stab,order,'blue');
        else
            stab_plot(om,Aout,stab,order,'blue','nolegends');
        end
        % plot special points
        SNidx = FRCom{k}.SNidx;
        SNfig = plot(om(SNidx),Aout(SNidx),thm.SN{:});
        set(get(get(SNfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        HBidx = FRCom{k}.HBidx;
        HBfig = plot(om(HBidx),Aout(HBidx),thm.HB{:});
        set(get(get(HBfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');         
    end    
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    zk = strcat('$||z_{',num2str(outdof(i)),'}||_{\infty}$');
    ylabel(zk,'Interpreter','latex');
    set(gca,'FontSize',14);
    grid on, axis tight; legend boxoff;
    
    % 3D plot
    figure; hold on
    for k=1:numLabs
        om = FRCom{k}.om;
        epsf = FRCom{k}.ep;
        Aout = FRCom{k}.Aout_frc(:,i);
        stab = FRCom{k}.st;
        if any(legDisp==k)
            stab_plot3(om,epsf,Aout,stab,order);
        else
            stab_plot3(om,epsf,Aout,stab,order,'nolegends');
        end
        % plot special points
        SNidx = FRCom{k}.SNidx;
        SNfig = plot3(om(SNidx),epsf(SNidx),Aout(SNidx),thm.SN{:});
        set(get(get(SNfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        HBidx = FRCom{k}.HBidx;
        HBfig = plot3(om(HBidx),epsf(HBidx),Aout(HBidx),thm.HB{:});
        set(get(get(HBfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');        
    end
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\epsilon$','interpreter','latex','FontSize',16);
    zk = strcat('$||z_{',num2str(outdof(i)),'}||_{\infty}$');
    zlabel(zk,'Interpreter','latex');
    set(gca,'FontSize',14);
    grid on, axis tight; legend boxoff;
    view([-1,-1,1]);
end
end
end


