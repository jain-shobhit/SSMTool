function bds = extract_FRC(obj, omega_range, varargin)
% EXTRACTFRC This function extracts the Frequency Response Curve in a given
% frequency span using continution package COCO
%
% bds = EXTRACT_FRC_FROMIC(obj, omega_range, varargin)
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
% See also: EXTRACT_FRC_FROMIC

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
omega0 = omega_range(1);
T0     = 2*pi/omega0*obj.periodsRatio;
tf     = obj.nCycles*T0;
outdof = obj.outdof;
epsilon = obj.system.Fext.epsilon;

if isempty(obj.system.Omega)
    obj.system.Omega = omega0;
end
assert(numel(obj.system.Omega)==1, 'coco run assumes single freq component');

switch obj.initialGuess
    case 'forward'
        %% initial solution by forward simulation
        % ode45 is used here. Integration option may be added in future
        x0_init = zeros(N,1);
        odefun = @(t,x) obj.ode_het(t,x,[omega0;epsilon]);
        [~, x0_po] = ode15s(odefun , [0 tf], x0_init);                % transient
        options = odeset('RelTol', 1e-9, 'AbsTol',1e-9);
        [t0, x0] = ode45(odefun, [0 T0], x0_po(end,:)', options);    % steady state
    case 'linear'
        %% initial solution by solving linear equation M\ddot{x}+C\dot{x}+Kx = F cos(O*t)
        if obj.system.order == 2
            mass = obj.system.M;
            damp = obj.system.C;
            stif = obj.system.K;
            fext = 2*obj.system.fext.data(1).f_n_k(1).coeffs;
            fext = epsilon*fext;
            kapa = obj.system.fext.data(1).kappa;
            if kapa<0
                kapa = -kapa;
            end
            qcom = (-(kapa*omega0)^2*mass+1i*kapa*omega0*damp+stif)\fext(:,1);
            t0   = linspace(0,T0,100);
            solx = qcom*exp(1i * kapa * omega0 * t0);
            solv = 1i * kapa * omega0 *solx;
            x0   = real([solx; solv])';
        else
            Fext = 2*obj.system.Fext.data(1).F_n_k(1).coeffs;
            Fext = epsilon*Fext;
            kapa = obj.system.Fext.data(1).kappa;
            if kapa<0
                kapa = -kapa;
            end            
            x0 = (1i*kapa*omega0*obj.system.B-obj.system.A)\Fext;
            t0   = linspace(0,T0,100);
            x0 = x0*exp(1i * kapa * omega0 * t0);
            x0 = real(x0)';
        end
end

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
bd0    = coco(prob, runid, [], cont_args{:});

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
% plot results
figure(gcf); hold on
if numoutdof>1
    for k=1:numoutdof
        subplot(numoutdof,1,k); hold on
        if ~epsflag
            plot(omega(stab), amp{k}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
            plot(omega(~stab), amp{k}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
        else
            plot(epsf(stab), amp{k}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
            plot(epsf(~stab), amp{k}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
        end
    end
else
    if ~epsflag
        plot(omega(stab), amp{1}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
        plot(omega(~stab), amp{1}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
    else
        plot(epsf(stab), amp{1}(stab), 'g*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-stable');
        plot(epsf(~stab), amp{1}(~stab), 'c*', 'MarkerSize', 6,'LineWidth',2,'DisplayName','COCO-unstable');
    end        
end
legend('show');
end
% figure;
% thm = struct( 'special', {{'FP', 'PD', 'TR'}});
% thm.FP = {'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'cyan', 'Marker', 'v', 'MarkerSize', 10, 'MarkerEdgeColor', ...
%   'cyan', 'MarkerFaceColor', 'white'};
% thm.PD = {'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'black', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
%   'black', 'MarkerFaceColor', 'white'};
% thm.TR = {'LineStyle', 'none', 'LineWidth', 2, ...
%   'Color', 'red', 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', ...
%   'red', 'MarkerFaceColor', 'white'};
% for i=1:numBP+1
%     if isempty(varargin)
%         coco_plot_bd(thm, runids{i}, 'omega', ampNames{1})
%     else
%         coco_plot_bd(thm, runids{i}, 'eps', ampNames{1})
%     end
% end
% grid on; box on;
% set(gca,'LineWidth',1.2);
% set(gca,'FontSize',14);
% if isempty(varargin)
% xlabel('$\Omega$','interpreter','latex','FontSize',16);
% else
% xlabel('$\epsilon$','interpreter','latex','FontSize',16);
% end
% ylabel(strcat('$||z_{',num2str(outdof(1)),'}||_{\infty}$'),'interpreter','latex','FontSize',16);
% title('FRC by coco(solid/dashed - stable/unstable)');
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
