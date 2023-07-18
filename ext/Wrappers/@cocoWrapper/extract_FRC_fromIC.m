function bds = extract_FRC_fromIC(obj, omega_range,IC,omegaIC, varargin)
% EXTRACT_FRC_FROMIC This function extracts the Frequency Response Curve in a given
% frequency span using continution package COCO and starting with a given
% initial condition. This enables continuation of specific branches of
% response
%
% bds = EXTRACT_FRC_FROMIC(obj, omega_range,IC,omegaIC, varargin)
%
% obj:      SSM class object
% omega_range:
%           range over which FRC is computed
% IC:       initial condition of periodic orbit
% omegaIC:  frequency of the initial periodic orbit
% varargin: flags to determine whether continuation is to be carried out in
%           epsilon or omega and if FRC is to be plotted
%
% bds:      bifurcation data struct that contains the information on the
%           continuation problems that are investigated
%
% See also: EXTRACT_FRC

% flag for parName and plotting
if numel(varargin)==2
        epsflag = true;
        PlotFRC = varargin{2};
elseif numel(varargin)==1
    if islogical(varargin{1})
        PlotFRC = varargin{1};
        epsflag = false;
    else
        PlotFRC = true;
        epsflag = true;
    end
else
    epsflag = false;
    PlotFRC = true;
end

%% initial setup
dir_name = obj.Options.dir_name;
N = obj.system.N;
n = obj.system.n;
omega0 = omegaIC;
T0     = 2*pi/omega0*obj.periodsRatio;
tf     = obj.nCycles*T0;
outdof = obj.outdof;
epsilon = obj.system.Fext.epsilon;

if isempty(obj.system.Omega)
    obj.system.Omega = omega0;
end
assert(numel(obj.system.Omega)==1, 'coco run assumes single freq component');

% Initial Guess
odefun = @(t,x) obj.ode_het(t,x,[omega0;epsilon]);

options = odeset('RelTol', 1e-9, 'AbsTol',1e-9);
[t0, x0] = ode45(odefun, [0 T0], IC.', options);    % steady state

%% continuation excitation frequency
% setup coco
prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
odefun = @(t,x,p) obj.ode_het(t,x,p);
funcs = {odefun};

coll_args = {funcs{:}, t0, x0, {'omega','eps'}, [omega0,epsilon]};   %#ok<CCAT>
prob = ode_isol2po(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
omData = struct();
omData.periodsRatio = obj.periodsRatio;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);

% track amplitude of outdof
ampdata.dof  = outdof;
ampdata.zdim = N;
numoutdof = numel(outdof);
ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
   ampNames{k} = strcat('amp',num2str(outdof(k)));
end
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
if ~epsflag
    cont_args = {1, {'omega' 'po.period' 'eps' ampNames{:}}, [omega_range(1) omega_range(end)]};
else
    eps_range = varargin{1};
    cont_args = {1, {'eps' 'omega' ampNames{:}}, [eps_range(1) eps_range(2)]};
end
runname = strcat('FRC',num2str(epsilon));

runid = coco_get_id(dir_name, runname);
fprintf('\n Run=''%s'': Continue primary family of periodic orbits.\n', ...
  runid);


%% 0 dim run to find Periodic Orbit with given Frequency
cont_args = {0, {'omega' 'po.period' 'eps' ampNames{:}}, [omega_range(1) omega_range(end)]};

bdinit    = coco(prob, runid, [], cont_args{:});

%% normal run
if ~epsflag
    cont_args = {1, {'omega' 'po.period' 'eps' ampNames{:}}, [omega_range(1) omega_range(end)]};
else
    eps_range = varargin{1};
    cont_args = {1, {'eps' 'omega' ampNames{:}}, [eps_range(1) eps_range(2)]};
end

EPlab = coco_bd_labs(bdinit, 'EP');
runname = strcat('FRC0',num2str(epsilon));
runid0 = coco_get_id(dir_name, runname);

prob = coco_prob();
prob = cocoSet(obj, prob);
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'ode', 'vectorized', true);
prob = ode_po2po(prob, '', runid, EPlab(1));
[data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
    'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
    'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
bd0 = coco(prob, runid0, [], cont_args{:});

% continuation along secondary branch if branch points are detected along
% the continuation above.

if obj.branchSwitch
    BPlab = coco_bd_labs(bd0, 'BP');
    numBP = numel(BPlab);
    bds = cell(numBP+1,1);
    runids = cell(numBP+1,1);
    bds{1} = bd0;
    runids{1} = runid;
    for k=1:numBP
        runidk = coco_get_id(runid, num2str(k));
        prob = coco_prob();
        prob = cocoSet(obj, prob);
        prob = coco_set(prob, 'ode', 'autonomous', false);
        prob = coco_set(prob, 'ode', 'vectorized', true);
        prob = ode_BP2po(prob, '', runid, BPlab(k));
        [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
        maps = data.coll_seg.maps;
        prob = coco_add_func(prob, 'OmegaT', @OmegaT, @OmegaT_du, omData, 'zero',...
            'uidx', [uidx(maps.T_idx), uidx(maps.p_idx(1))]);
        prob = coco_add_func(prob, 'amp', @amplitude, ampdata, 'regular', ampNames,...
            'uidx', uidx(maps.xbp_idx), 'remesh', @amplitude_remesh);
        fprintf('\n Run=''%s'': Continue equilibria along secondary branch.\n', ...
                  runidk);
        bdk = coco(prob, runidk, [], cont_args{:});
        bds{k+1} = bdk;
        runids{k+1} = runidk;
    end
else
    numBP = 0;
    bds = {bd0};
    runids = {runid};
end


% result visualization
% extract results
omega = [];
epsf  = [];
stab = logical([]);
amp = cell(numoutdof,1);
for i=1:numBP+1
    bd = bds{i};
    omegai = coco_bd_col(bd, 'omega');
    epsfi = coco_bd_col(bd, 'eps');
    stabi = coco_bd_col(bd, 'eigs')';
    stabi = abs(stabi);
    stabi = all(stabi<1, 2);
    omega = [omega, omegai];
    epsf  = [epsf, epsfi];
    stab  = [stab; stabi];
    for j=1:numoutdof
        ampj   = coco_bd_col(bd, ampNames{j});
        amp{j} = [amp{j}, ampj];
    end
end

if PlotFRC
    
    ustabc = [0.4660    0.6740    0.1880];
    stabc = [0.3010    0.7450    0.9330];
% plot results
figure(gcf); hold on
if numoutdof>1
    for k=1:numoutdof
        subplot(numoutdof,1,k); hold on
        if ~epsflag
            plot(omega(stab), amp{k}(stab), '*','color',stabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
            plot(omega(~stab), amp{k}(~stab), '*','color',ustabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
        else
            plot(epsf(stab), amp{k}(stab), '*','color',stabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
            plot(epsf(~stab), amp{k}(~stab), '*','color',ustabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
        end
    end
else
    if ~epsflag
        plot(omega(stab), amp{1}(stab), '*','color',stabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
        plot(omega(~stab), amp{1}(~stab), '*','color',ustabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
    else
        plot(epsf(stab), amp{1}(stab), '*','color',stabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
        plot(epsf(~stab), amp{1}(~stab), '*','color',ustabc, 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
    end
end
legend('show');
end
end



function [data, y] = OmegaT(prob, data, u) %#ok<INUSL>
y = u(1)*u(2)-2*pi*data.periodsRatio; % omega*T = 2*pi
end

function [data, J] = OmegaT_du(prob, data, u) %#ok<INUSL>
J = [u(2) u(1)];
end

function [data, y] = amplitude(prob, data, u) %#ok<INUSL>

xbp = reshape(u, data.zdim, []);
y = xbp(data.dof,:);
y = max(abs(y),[],2);

end

function [prob, status, xtr] = amplitude_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[colldata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = colldata.coll_seg.maps;
xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
status = 'success';

end



function y = ode_het_dfdx(t, x, p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = -x2-p1.*x1+cos(t);

end

function y = ode_het_dfdp(t, x, p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = -x2-p1.*x1+cos(t);

end

function y = ode_het_dfdt(t, x, p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = -x2-p1.*x1+cos(t);

end
