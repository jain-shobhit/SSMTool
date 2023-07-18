function [SD] = extract_Stability_Diagram(obj,omRange,epsRange, parName,p0,varargin)
% EXTRACT_STABILITY_DIAGRAM This function extracts the stability diagram
% using the full dynamical system by continuing either PD or SN
% bifurcations. These families of bifurcations form the boundary regions of
% stable and unstable domains of the trivial response of parametrically
% excited systems. The resulting graphs are also known as stability diagram
% or Ince Strutt diagrams.
%
% FRC = EXTRACT_STABILITY_DIAGRAM(OBJ,OMRANGE,EPSRANGE,PARNAME,PO,VARARGIN)
%
% omRange :     range of forcing frequency over which stability is evaluated
% epsRange:     range of forcing amplitude over which stability is evaluated
% parName :     'amp' / 'freq' determines which parameter is released initially
%               when looking for bifurcations. 
% p0      :     [Omega_i, Epsilon_i] Initial parameter values for
%               continuation
% varargin:     {bifType , PlotSD} bifType = 'SN', 'PD' whether stability diagram 
%               is obtained by continuingsaddle node or period doubling bifurcations
%               PlotSD = true/false can be used to suppress plot
%
% SD:           Stability Diagram data struct
%
% See also: SSM/EXTRACT_STABILITY_DIAGRAM

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
%% ODE Function

odefun = @(t,z,p) ode_het(obj, t, z, p);
funcs = {odefun};

%% Initial Condition
switch bifType
    case 'PD'
        periodsRatio = 1;
    case 'SN'
        periodsRatio = 2;
end

T = periodsRatio*2*pi/p0(1); % Forcing period of initial condition

t0 = linspace(0,T,101);
x0 = zeros(101,obj.system.N);

%% Build continuation problem structure

prob = coco_prob();
prob = cocoSet(obj,prob);
prob = coco_set(prob, 'ode', 'autonomous', false,'vectorized', true);

coll_args = [funcs, {t0, x0, {'om' 'eps' }, p0}];
prob = ode_isol2po(prob, '', coll_args{:});

% Fix period of response to forcing frequency to detect PD bifurcation
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;

omData = struct();
omData.periodsRatio = periodsRatio;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

%% Continuation
switch parName
    case 'freq'
        
        isomega = true;
        cont_args = { 1, { 'om', 'po.period' }, omRange };
        parRange = epsRange; %Range of 2nd continuation parameter
    case 'amp'
        isomega = false;
        cont_args = {1, {'eps','om'}, epsRange };
        parRange = omRange; %Range of 2nd continuation parameter
    otherwise
        error('Continuation parameter should be freq or amp');
end

runid = 'full_detect_bif';
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', runid);

bd = coco(prob,runid,[],cont_args{:});


%% Continue PD branches
switch bifType
    case 'PD'
        labs = coco_bd_labs(bd, 'PD');
    case 'SN'
        labs = coco_bd_labs(bd, 'SN');
end

if obj.branchSwitch
    if isempty(labs)
        sprintf('No bifurcations were detected');
        SD = [];
        return
    else
        Epsilon = [];
        Omega   = [];
        run_idx = 1;
        for lab = labs
            runid_bif = sprintf('full_family_bif_%d',run_idx ); % Id of this run along family of PD bifurcations
            fprintf('\n Run=''%s'': Continue bifurcations from point %d in run ''%s''.\n',runid_bif, lab, runid);
            
            [omega_i,epsilon_i] = continue_bif(obj,parRange,lab, runid,runid_bif,isomega,omData,bifType);
            Omega   = [Omega,omega_i];
            Epsilon = [Epsilon,epsilon_i];
            run_idx = run_idx+1;
        end
        
        
        [SD.Omega,idx] = sort(Omega);
        SD.Epsilon = Epsilon(idx);
    end
    
    if PlotSD
        plot_SD(SD, epsRange)
    end
    
else
    SD = [];
end
end

function[omega,epsilon] = continue_bif(obj,parRange,lab, labid,runid,isomega,omData,bifType)
% For continuation of families of bifurcations in omega and epsilon

prob = coco_prob();
prob = cocoSet(obj,prob); % causes issues when used with PD

switch bifType
    case 'PD'
        prob = ode_PD2PD(prob, '', labid, lab);
    case 'SN'
        prob = ode_SN2SN(prob, '', labid, lab);
end

[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
            'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);


if isomega
    cont_args = {1, {'eps' 'om' 'po.period' }, parRange};
else
    cont_args = {1, {'om' 'po.period' 'eps' }, parRange};
end

bd  = coco(prob, runid, [], cont_args{:});

% read out omega and epsilon values along period doubling bifurcation
omega = coco_bd_col(bd, 'om');
epsilon = coco_bd_col(bd, 'eps');
end


function plot_SD(SD, epsRange)
figure(gcf); hold on;box on
ymax = max(SD.Epsilon) * 1.1;

%% Full model SD
plot(SD.Omega, SD.Epsilon,'ok', 'MarkerSize', 5, 'DisplayName', 'Full System')

%% Shading of the Figure
%{
a = area(SD.Omega,[SD.Epsilon ; ones(1,numel(SD.Omega))*ymax - SD.Epsilon].' ,'FaceAlpha',.3);

a(1).FaceColor = [0.9 0.9 0.9];
a(1).LineWidth = 0.001;
a(1).DisplayName = 'Stable';

a(2).FaceColor = [0.5 0.5 0.5];
a(2).LineWidth = 0.001;
a(2).DisplayName = 'Unstable';
%}
%% Plot Specs
add_labels('$\Omega$','$\epsilon$')
add_legends();

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%title('Stability of trivial fixed point')

% Y axis limit
ylim([epsRange(1),ymax]);

end

function add_legends()
legend('show');
legend boxon;
end

function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');

grid on, axis tight; legend boxoff;
end

function [data, y] = OmegaT(prob, data, u) %#ok<INUSL>
y = u(1)*u(2)-2*pi*data.periodsRatio; % omega*T = 2*pi
end

function [data, J] = OmegaT_du(prob, data, u) %#ok<INUSL>
J = [u(2) u(1)];
end


