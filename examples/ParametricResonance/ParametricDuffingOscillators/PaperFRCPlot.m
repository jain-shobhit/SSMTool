function PaperFRCPlot(FRC,bd,order,outdof,varargin)
% FRC data from SSM
% bd data from coco
numOutdof = numel(outdof);


%% Prepare SSM data
Par   = [FRC.Omega];
Aout  = [FRC.Aout];
stab_ssm  = [FRC.stability];

%numPts    = numel(stab_ssm);
%Aout = reshape(Aout, [numOutdof, numPts]);
%Aout = Aout';

%% Prepare COCO data
nSample = 1;

bd = bd{1};
% extract results
omega = [];
stab  = logical([]);
amp   = cell(numOutdof,1);
omegai = coco_bd_col(bd, 'omega');
stabi  = coco_bd_col(bd, 'eigs')';
stabi  = abs(stabi);
stabi  = all(stabi<1, 2);
omega  = [omega, omegai];
stab   = [stab; stabi];

ampNames = cell(1, numel(numOutdof));
for k = 1:numel(outdof)
    ampNames{k} = strcat('amp',num2str(outdof(k)));
end

for j=1:numOutdof
    ampj   = coco_bd_col(bd, ampNames{j});
    amp{j} = [amp{j}, ampj];
end

if ~isempty(varargin)
    % Second Coco bd input
    bd2 = varargin{1};
    bd2 = bd2{1};
    omega2 = [];
    stab2  = logical([]);
    amp2   = cell(numOutdof,1);
    omegai2 = coco_bd_col(bd2, 'omega');
    stabi2  = coco_bd_col(bd2, 'eigs')';
    stabi2  = abs(stabi2);
    stabi2  = all(stabi2<1, 2);
    omega2  = [omega2, omegai2];
    stab2   = [stab2; stabi2];
    
    ampNames2 = cell(1, numel(numOutdof));
    for k = 1:numel(outdof)
        ampNames2{k} = strcat('amp',num2str(outdof(k)));
    end
    
    for j=1:numOutdof
        ampj2   = coco_bd_col(bd2, ampNames2{j});
        amp2{j} = [amp2{j}, ampj2];
    end
    
end

%% Plotting parameters

%% Single Plot Setting

%axis_size = 20; % Set legend and Axis size
%xwidth = 500; % window size
%ywidth = 400;

%ssm_marker = 10*1.4;
%coco_marker = 6*1.4;


%% Double Plot Setting
% {
axis_size = 20; % Set legend and Axis size
xwidth = 600; % window size
ywidth = 500;

ssm_marker = 2*1.4;
coco_marker = 6*1.4;
%}

ustabc1 = [0.9290    0.6940    0.1250];
stabc1 = [0.4940    0.1840    0.5560];
ustabc2 = [0.4660    0.6740    0.1880];
stabc2 = [0.3010    0.7450    0.9330];
%% Start Plotting
for k=1:numOutdof
    set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
    figure(); hold on

    % Plot SSM Data    
    plot(Par(stab_ssm),Aout(stab_ssm),'bo','MarkerSize',ssm_marker);
    plot(Par(~stab_ssm),Aout(~stab_ssm),'or','MarkerSize',ssm_marker);
    add_labels('$\Omega$',strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$'))
    
    % Plot Coco Data
    om = omega(~stab); am = amp{k}(~stab);
    plot(om(1:nSample:end), am(1:nSample:end), '*', 'MarkerSize',coco_marker,'color',ustabc1);
    
    om = omega(stab); am = amp{k}(stab);
    plot(om(1:nSample:end), am(1:nSample:end), '*', 'MarkerSize',coco_marker,'color',stabc1);

    
    % Plot second set of Coco data
    if ~isempty(varargin)
        om2 = omega2(~stab2); am2 = amp2{k}(~stab2);
        plot(om2(1:nSample:end), am2(1:nSample:end), '*', 'MarkerSize',coco_marker,'color',ustabc2);
        % Plot Coco Data
        om2 = omega2(stab2); am2 = amp2{k}(stab2);
        plot(om2(1:nSample:end), am2(1:nSample:end), '*', 'MarkerSize',coco_marker,'color',stabc2);
        
        %'color', [0.9294    0.6941    0.1255],
    end
    

    
    % Add legend
    if isempty(varargin)
        add_legends(stab_ssm,order, get_cocolabel(stab))
    else
            add_legends(stab_ssm,order, get_cocolabel(stab,stab2))

    end
    
        
    %% Set axis limits by SSM axis
    xmin = min(Par);
    xmax = max(Par);
    %ymin = min(Aout(:,k));
    %ymax = max(Aout(:,k))*1.1;
    
    xlim([xmin,xmax]);
    %ylim([ymin,ymax]);
    
    %% Parameters for figure size, label size, legend size
    set(gca,'fontsize',axis_size)
    
    %% Size of window
    set(gcf, 'Position', [100 100 xwidth ywidth])
end

end

function cclabel = get_cocolabel(stab,varargin)
cclabel = {'a'};


if all(stab)
    cclabel{end+1} = 'collocation-stable' ;
elseif ~any(stab)
    cclabel{end+1} = 'collocation-unstable';
else
    cclabel{end+1} = 'collocation-unstable';
    cclabel{end+1} = 'collocation-stable' ;
    
end

if ~isempty(varargin)
    stab2 = varargin{1};
    if all(stab2)
        cclabel{end+1} = 'collocation-stable' ;
    elseif ~any(stab2)
        cclabel{end+1} = 'collocation-unstable';
    else
        cclabel{end+1} = 'collocation-unstable';
        cclabel{end+1} = 'collocation-stable' ;
        
    end
    

end

cclabel = cclabel(2:end);

end

function add_legends(stab,order,cclabel)

if all(stab)
    legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),cclabel{:},...
        'Interpreter','latex','Location','best');
elseif ~any(stab)
    legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),cclabel{:},...
        'Interpreter','latex','Location','best');
else
    legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
        strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),cclabel{:},...
        'Interpreter','latex','Location','best');
end

end

function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');
set(gca,'FontSize',14);
grid on, axis tight; legend boxoff;
end