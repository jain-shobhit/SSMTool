function PBplotSweep(FRCom,bds,outdof,order,varargin)
figure()

%
%bds = vertcat(bds{:})

ssmwd = 1.5;
ccsz = 3*1.4;

%% results extraction
nSample = 2;
numLabstot = numel(bds);
FRCombd = cell(numLabstot,1);

ampNames = cell(1, numel(outdof));
for i = 1:numel(outdof)
    ampNames{i} = strcat('amp',num2str(outdof(i)));
end

for k=1:numLabstot
    
    bdk = bds{k}; 
    
    stabi = coco_bd_col(bdk, 'eigs')';
    stabi = abs(stabi);
    
    FRCombd{k}.om = coco_bd_col(bdk, 'omega');
    FRCombd{k}.ep = coco_bd_col(bdk, 'eps');
    FRCombd{k}.st = all(stabi<1, 2);
    
    
    
    for j=1:numel(outdof)
        ampj   = coco_bd_col(bdk, ampNames{j});
        FRCombd{k}.Aout_frc(:,j)= ampj.';
    end
    
    if  any(FRCombd{k}.Aout_frc) % sparsify
        FRCombd{k}.om = FRCombd{k}.om(1:nSample:end);
        FRCombd{k}.ep = FRCombd{k}.ep(1:nSample:end);
        FRCombd{k}.st = FRCombd{k}.st(1:nSample:end) ;
        FRCombd{k}.Aout_frc = FRCombd{k}.Aout_frc(1:nSample:end);
    end
    
end



sweeps_plotcc(numLabstot,FRCombd,outdof,ccsz);
%% SSMplot
numLabstot_ssm = numel(FRCom);
sweeps_plot(numLabstot_ssm,FRCom,outdof,order,ssmwd);
%% Plot Settings

axis_size = 20; % Set legend and Axis size
xwidth = 600; % window size
ywidth = 500;

%% Parameters for figure size, label size, legend size
set(gca,'fontsize',axis_size)

%% Size of window
set(gcf, 'Position', [0 0 xwidth ywidth])

y_tick = [];
for k = 1:numel(FRCom)
    epsf = FRCom{k}.ep(1);
    y_tick = [y_tick;epsf];
end
y_tick = unique(y_tick);
yticks(y_tick)
ylim([0.18,0.62])
grid on, axis tight;legend boxoff;
view([-1,-1,1]);
end


function sweeps_plotcc(numLabs,FRCom,outdof,ccsz)

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

ndof = numel(outdof);
for i=1:ndof
    % 3D plot
    figure(gcf); hold on
    for k=1:numLabs
        om = FRCom{k}.om;
        epsf = FRCom{k}.ep;
        Aout = FRCom{k}.Aout_frc(:,i);
        stab = FRCom{k}.st;
        if any(legDisp==k)
            stab_plot3cc(om,epsf,Aout,stab,ccsz);
        else
            stab_plot3cc(om,epsf,Aout,stab,ccsz,'nolegends');
        end
    end
    
end

end

function stab_plot3cc(x,y,z,S,ccsz,varargin)
% This function is adapted from stab_plot in coco_plot_bd in coco

II  = 1;
EI  = numel(x);
I   = II;
figs = [];
stab = [];
ST = cell(2,1);
ST{1} = {'c*','MarkerSize', ccsz};
ST{2} = {'g*','MarkerSize', ccsz};
%ST{1} = {'b--','LineWidth',1.5};
%ST{2} = {'b-','LineWidth',1.5};
while true
    if II>=EI; break; end
    [I, I0] = next_index(EI, S, I, II);
    fig = plot3(x(II:I0), y(II:I0), z(II:I0), ST{S(II)+1}{:});
    figs = [figs fig];
    stab = [stab S(II)+1];
    II = I0;
end

if isempty(varargin)
    Leg = cell(2,1);
    Leg{1} = strcat('Collocation - unstable');
    Leg{2} = strcat('Collocation - stable');
    numSegs = numel(figs);
    if numSegs==1
        % add legend
        legend(Leg{stab(1)},'Interpreter','latex');
    else
        legend(Leg{stab(1)},Leg{stab(2)},'Interpreter','latex');
        if numSegs>2
            % remove repetitive legends from option
            for j=3:numSegs
                set(get(get(figs(j),'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
            end
        end
    end
else
    numSegs = numel(figs);
    for j=1:numSegs
        set(get(get(figs(j),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
end
end

function sweeps_plot(numLabs,FRCom,outdof,order,ssmsz)
% Plot FRC in normal coordinates
thm = struct( 'special', {{'SN', 'HB' , 'PD'}});
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
    'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
    'black', 'MarkerFaceColor', 'white'};
thm.PD = {'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'green', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
    'black', 'MarkerFaceColor', 'white'};



% Plot FRC in system coordinates
% setup legend
legTypes = zeros(numLabs,1); % each entry can be 1/2/3 - all stab/all unstab/mix
for k=1:numLabs
    stab = FRCom{k}.stability;
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

ndof = numel(outdof);
for i=1:ndof
    
    hold on
    for k=1:numLabs
        om = FRCom{k}.om;
        epsf = FRCom{k}.ep;
        Aout = FRCom{k}.Aout(i,:);
        stab = FRCom{k}.stability;
        if any(legDisp==k)
            stab_plot3(om,epsf,Aout,stab,order,ssmsz);
        else
            stab_plot3(om,epsf,Aout,stab,order,ssmsz,'nolegends');
        end
        % plot special points
        SNidx = FRCom{k}.SNidx;
        SNfig = plot3(om(SNidx),epsf(SNidx),Aout(SNidx),thm.SN{:});
        set(get(get(SNfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        %HBidx = FRCom{k}.HBidx;
        %HBfig = plot3(om(HBidx),epsf(HBidx),Aout(HBidx),thm.HB{:});
        %set(get(get(HBfig,'Annotation'),'LegendInformation'),...
        %    'IconDisplayStyle','off');
    end
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\epsilon$','interpreter','latex','FontSize',16);
    zk = strcat('$||z_{',num2str(outdof(i)),'}||_{\infty}$');
    zlabel(zk,'Interpreter','latex');
    
    grid on, axis tight; legend boxoff;
    view([-1,-1,1]);
end

end

function stab_plot3(x,y,z,S,order,ssmsz,varargin)
% This function is adapted from stab_plot in coco_plot_bd in coco

II  = 1;
EI  = numel(x);
I   = II;
figs = [];
stab = [];
ST = cell(2,1);
%ST{1} = {'ro','LineWidth',ssmsz};
%ST{2} = {'bo','LineWidth',ssmsz};
ST{1} = {'b--','LineWidth',ssmsz};
ST{2} = {'b-','LineWidth',ssmsz};
while true
    if II>=EI; break; end
    [I, I0] = next_index(EI, S, I, II);
    fig = plot3(x(II:I0), y(II:I0), z(II:I0), ST{S(II)+1}{:});
    figs = [figs fig];
    stab = [stab S(II)+1];
    II = I0;
end

if isempty(varargin)
    Leg = cell(4,1);
    Leg{1} = strcat('Collocation -  stable');
    Leg{2} = strcat('Collocation - unstable');
    Leg{3} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable');
    Leg{4} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable');
    numSegs = numel(figs);
    if numSegs==1
        % add legend
        legend(Leg{stab(1)},'Interpreter','latex');
    else
        legend(Leg{:}, 'Interpreter','latex');
        if numSegs>4
            % remove repetitive legends from option
            for j=5:numSegs
                set(get(get(figs(j),'Annotation'),'LegendInformation'),...
                    'IconDisplayStyle','off');
            end
        end
    end
else
    numSegs = numel(figs);
    for j=1:numSegs
        set(get(get(figs(j),'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
    end
end
end

function [I, I0] = next_index(EI, S, I, II)
while (I<EI && S(II)==S(I)); I=I+1; end
I0 = I;
end

%{
%% Double Plot Setting
% {
axis_size = 20; % Set legend and Axis size
xwidth = 600; % window size
ywidth = 500;

ssmwd = 1.5;
cocosz = 7 *1.4;
% }

ustabc1 = [0.9290    0.6940    0.1250];
stabc1 = [0.4940    0.1840    0.5560];
ustabc2 = [0.4660    0.6740    0.1880];
stabc2 = [0.3010    0.7450    0.9330];
%% Prepare SSMdata
sSample = 1;
numoutdof = numel(outdof);



%% Prepare COCOdata
nSample = 1;
inSample = 1;
ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
    ampNames{k} = strcat('amp',num2str(outdof(k)));
end

omegacc = [];
epsfcc  = [];
stabcc = logical([]);
amp = cell(numoutdof,1);

for i=1:numel(bds)
    bd = bds{i};
    bd = bd{1};
    omegai = coco_bd_col(bd, 'omega');
    epsfi = coco_bd_col(bd, 'eps');
    stabi = coco_bd_col(bd, 'eigs')';
    stabi = abs(stabi);
    stabi = all(stabi<1, 2);
    omegacc = [omegacc, omegai];
    epsfcc  = [epsfcc, epsfi];
    stabcc  = [stabcc; stabi];
    for j=1:numoutdof
        ampj   = coco_bd_col(bd, ampNames{j});
        amp{j} = [amp{j}, ampj];
    end
end

if ~isempty(varargin)
    iomegacc = [];
    iepsfcc  = [];
    istabcc = logical([]);
    iamp = cell(numoutdof,1);

    
    ibd = varargin{1};
    iomegai = coco_bd_col(ibd, 'omega');
    iepsfi = coco_bd_col(ibd, 'eps');
    istabi = coco_bd_col(ibd, 'eigs')';
    istabi = abs(istabi);
    istabi = all(istabi<1, 2);
    iomegacc = [iomegacc, iomegai];
    iepsfcc  = [iepsfcc, iepsfi];
    istabcc  = [istabcc; istabi];
    iamp{1} = [ibd{2:end,18}];
end

Leg = {'anfang'};
for i=1:numoutdof
    figure()
    hold on
    
    
        
    %% Coco Plot
%{
    om = omegacc(stabcc); am = amp{k}(stabcc); eps = epsfcc(stabcc);
    plot3(om(1:nSample:end),  eps(1:nSample:end),am(1:nSample:end), 'g*','MarkerSize', cocosz);%,...
        %'color',[0 0 0] );
    om = omegacc(~stabcc); am = amp{k}(~stabcc);eps = epsfcc(~stabcc);
    plot3(om(1:nSample:end), eps(1:nSample:end),am(1:nSample:end), 'c*', 'MarkerSize', cocosz);%, ...
        %'color',[0.9294    0.6941    0.1255] );
%}
%{
    om = omegacc(~stabcc); am = amp{k}(~stabcc);eps = epsfcc(~stabcc);
    plot3(om(1:nSample:end), eps(1:nSample:end),am(1:nSample:end), '*', 'MarkerSize', cocosz,'color',ustabc1)%, 'color',[0.9294    0.6941    0.1255] );
    om = omegacc(stabcc); am = amp{k}(stabcc); eps = epsfcc(stabcc);
    plot3(om(1:nSample:end),  eps(1:nSample:end),am(1:nSample:end), '*','MarkerSize', cocosz,'color',stabc1)%,'color',[0 0 0] );

%}
     Leg{end+1} = strcat('collocation -unstable');
    Leg{end+1} = strcat('collocation-stable');
   

    if ~isempty(varargin)
        iom = iomegacc(~istabcc); iam = iamp{k}(~istabcc);ieps = iepsfcc(~istabcc);
        plot3(iom(1:inSample:end), ieps(1:inSample:end),iam(1:inSample:end), '*', 'MarkerSize', cocosz,'color',ustabc2)%, 'color',[0.9294    0.6941    0.1255] );
  
        iom = iomegacc(istabcc); iam = iamp{k}(istabcc); ieps = iepsfcc(istabcc);
        plot3(iom(1:inSample:end),  ieps(1:inSample:end),iam(1:inSample:end), '*','MarkerSize', cocosz,'color',stabc2)%,'color',[0 0 0] );
        Leg{end+1} = strcat('collocation -unstable');
        Leg{end+1} = strcat('collocation-stable');

    end
        
    %% SSM Plot
    
    %% Annahme nur 1 outdof
    
    for j = 1:numel(FRCom)
        epsfssm = [FRCom{j}.ep];
        stab  = [FRCom{j}.st];
        Aouti = [FRCom{j}.Aout_frc];
        omssm = [FRCom{j}.om];
        
        
        oms = omssm(stab); eps = epsfssm(stab); am = Aouti(stab);
        plot3(oms(1:sSample:end),  eps(1:sSample:end),am(1:sSample:end),'b-','LineWidth',ssmwd);
        oms  = omssm(~stab); eps = epsfssm(~stab); am = Aouti(~stab);
        plot3(oms(1:sSample:end),  eps(1:sSample:end),am(1:sSample:end),'b--','LineWidth',ssmwd);
        
    end
    Leg{end+1} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable');
    Leg{end+1} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable');
    
    %% Labels and Legend
    legend(Leg{2:end},...
        'Interpreter','latex','Location','best');
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\epsilon$','interpreter','latex','FontSize',16);
    zk = strcat('$||z_{',num2str(outdof(i)),'}||_{\infty}$');
    zlabel(zk,'Interpreter','latex');
    %% Parameters for figure size, label size, legend size
    set(gca,'fontsize',axis_size)
    
    %% Size of window
    set(gcf, 'Position', [0 0 xwidth ywidth])
    
    grid on, axis tight;legend boxoff;
    view([-1,-1,1]);
end
end
%}
