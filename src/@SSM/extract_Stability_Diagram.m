function SD = extract_Stability_Diagram(obj,resModes, order, omRange, epsRange, parName, p0,varargin)
% EXTRACT_STABILITY_DIAGRAM This function extracts the
% stability diagram of the trivial fixed point of parametrically excited
% systems (also known as InceStrutt diagram) in  a range of forcing frequency
% around the principal subharmonic 2:1 resonance using the reduced order model computed
% by the SSM. An appropriate SSM is constructed based on the resonant
% spectrum of system. Using continuation periodic orbits are detected and
% period doubling / saddle node bifurcations which appear at the boundaries of stable
% regions around these resonances are tracked and continued.
%
% IC = EXTRACE_STABILITY_DIAGRAM(OBJ, RESMODES, ORDER, OMRANGE, EPSRANGE, PARNAME, P0, VARARGIN)
%
%
% resModes:     Modes of the Eigenspace over which SSM is computed
% order   :     order of SSM epansion
% omRange :     range of forcing frequency over which stability diagram is
%               computed
% epsRange:     range of forcing amplitude over which stability diagram is
%               computed
% parName :     'amp' / 'freq' determines which parameter is released initially
%               when looking for bifurcations. 
% p0      :     Initial parameters for continuation
%               p0(1) contains forcing frequency, has to be in omRange
%               p0(2) contains epsilon, has to be in epsRange
% varargin:     { ['PD' or 'SN'], PlotSD}, whether the resonance tongue is to
%               obtained using PD or SN bifurcations, 
%               option to suppress plot of SD
%
% SD:           Stability Diagram data struct
%
% See also: COCO/EXTRACT_STABILITY_DIAGRAM

startStab = tic;
%% Parse input parameters

% flag for biftype and plotting 
if numel(varargin)==2
        SNflag = strcmp(varargin{1},'SN');
        PlotSD = varargin{2};
elseif numel(varargin)==1
    if islogical(varargin{1})
        PlotSD = varargin{1};
        SNflag = false;
    else
        PlotSD = true;
        SNflag = strcmp(varargin{1},'SN');
    end
else
    SNflag = false;
    PlotSD = true;
end

if SNflag
    bifType = 'SN';
else
    bifType = 'PD';
end
%% Compute Autonomous SSM
obj.choose_E(resModes);

% Initialise external forcing frequency
obj.System.Omega = p0(1);

% Compute SSM
[W, R] = obj.compute_whisker(order);

%% Reduced dynamics ODEfunction
% Data for odefunction

if obj.FRCOptions.omDepNonAuto % Assume strong dependence of coefficients on Omega
    
    fdata   = struct('order', order,'R',R,'W',W);
    odefun    = @(t,x,p) ode_2DSSM_cartesian(t,x,p,fdata,obj);
    odefun_dx = @(t,x,p) ode_2DSSM_cartesian_DFDX(t,x,p,fdata,obj);
    odefun_dp = @(t,x,p) ode_2DSSM_cartesian_DFDP(t,x,p,fdata,obj);
else
    
    % Compute Reduced dynamics and sensitivity coefficients
    obj.System.Omega = obj.FRCOptions.omDepNonAutoVal;
    [X, S] = obj.compute_perturbed_whisker(order-1,W,R);
    [~,DS] = obj.compute_sensitivity_coefficients(order-1,W,R,X);
    
    
    fdata   = struct('order', order,'R',R,'S',S,'DS',DS);
    odefun    = @(t,x,p) ode_2DSSM_cartesian_fixROM(t,x,p,fdata);
    odefun_dx = @(t,x,p) ode_2DSSM_cartesian_DFDX_fixROM(t,x,p,fdata);
    odefun_dp = @(t,x,p) ode_2DSSM_cartesian_DFDP_fixROM(t,x,p,fdata,obj.Options.contribNonAuto);
end

funcs  = {odefun,odefun_dx,odefun_dp};

%% Set up continuation problem
% Initial conditions

switch bifType
    case 'PD'
        periodsRatio = 1;
    case 'SN'
        periodsRatio = 2;
end

T = periodsRatio *2*pi/p0(1); % Forcing period of initial condition


t0 = linspace(0,T,101);
x0 = zeros(101,2);

% Build Problem
prob = coco_prob();
prob = cocoSet(obj.contOptions, prob); %set default options
prob = coco_set(prob, 'ode', 'autonomous', false,'vectorized', false);

coll_args = [funcs, {t0, x0, {'om' 'eps' }, p0}];
prob = ode_isol2po(prob, '', coll_args{:});

% Fix period of response to forcing frequency to detect PD bifurcation
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
omData = struct();
omData.periodsRatio = periodsRatio;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

%% Continuation arguments
switch parName
    case 'freq'
        
        isomega = true;
        cont_args = { 1, { 'om', 'po.period' }, omRange };
        parRange = epsRange; %Range of 2nd continuation parameter
    case 'amp'
        isomega = false;
        cont_args = {1, {'eps', 'po.period'}, epsRange };
        parRange = omRange; %Range of 2nd continuation parameter
    otherwise
        error('Continuation parameter should be freq or amp');
end

%% Continuation
runid = 'ROM_detect_bif';
bd = coco(prob,runid,[],cont_args{:});


%% Continue PD branches

switch bifType
    case 'PD'
        labs = coco_bd_labs(bd, 'PD');
    case 'SN'
        labs = coco_bd_labs(bd, 'SN');
end
if obj.FRCOptions.branchSwitch
    if isempty(labs)
        sprintf('No bifurcations indicating a resonance tongue were detected');
        SD = [];
        return
    else
        epsilon = [];
        omega   = [];
        run_idx = 1;
        for lab = labs
            
            labid = sprintf('ROM_family_bif%d',run_idx ); % Id of this run along family of PD bifurcations
            fprintf('\n Run=''%s'': Continue bifurcations from point %d in run ''%s''.\n',labid, lab, runid);
            
            [omega_i,epsilon_i] = stab_diag(obj,parRange,lab,runid,labid,isomega,bifType, omData);
            omega   = [omega,omega_i];
            epsilon = [epsilon,epsilon_i];
            run_idx = run_idx+1;
        end
        SD.omega = omega;
        SD.epsilon = epsilon;
    end
    
    if PlotSD 
        plot_Stability_Diagram(omega,epsilon,order);
    end
    
    SD.order = order;
    totalComputationTime = toc(startStab);
    disp(['Total time spent on Stability Diagram computation = ' datestr(datenum(0,0,0,0,0,totalComputationTime),'HH:MM:SS')])
else
    SD = [];    
end

end



function[omega,epsilon] = stab_diag(obj,parRange,lab, runid,labid,isomega, bifType, omData)
% For continuation of PD bifurcation in omega and epsilon

% Build continuation problem
prob = coco_prob();
prob = cocoSet(obj.contOptions, prob); %set default options
switch bifType
    case 'PD'
        prob = ode_PD2PD(prob, '', runid, lab);
    case 'SN'
        prob = ode_SN2SN(prob, '', runid, lab);
end

% Fix period of response
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
            'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);


if isomega
    cont_args = {1, {'eps' 'om' 'po.period' }, parRange};
else
    cont_args = {1, {'om' 'po.period' 'eps' }, parRange};
end
bd_lab  = coco(prob, labid, [], cont_args{:});

% read out omega and epsilon values along period doubling bifurcation

omega = coco_bd_col(bd_lab, 'om');
epsilon = coco_bd_col(bd_lab, 'eps');
end




function [] = plot_Stability_Diagram(omega,epsilon,order)
figSD = gcf;
% Plot in terms of omega and epsilon
figure(figSD); hold on; grid off;
plot(omega,epsilon,'-','LineWidth',2,'color',[0.4940    0.1840    0.5560],'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(order),')$$'))
add_labels('$\Omega$','$\epsilon$')
add_legends(order);
%title('Stability of trivial fixed point')
%legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$'),...
%        'Interpreter','latex','Location','best');
%legend boxon

end


function add_legends(order)
hl = legend('show','Location','best');
set(hl, 'Interpreter','latex')
legend boxon;
end


function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');
set(gca,'FontSize',14);
grid on, axis tight; legend boxoff;
end

function [data, y] = OmegaT(prob, data, u) %#ok<INUSL>
y = u(1)*u(2)-2*pi*data.periodsRatio; % omega*T = 2*pi
end

function [data, J] = OmegaT_du(prob, data, u) %#ok<INUSL>
J = [u(2) u(1)];
end


