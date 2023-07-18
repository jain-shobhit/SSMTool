function Sweepplot_for_paper_lvlset(FRCom,bds,outdof,order,mus,varargin)
%% Double Plot Setting
% {
axis_size = 20; % Set legend and Axis size
xwidth = 600; % window size
ywidth = 500;

ssmsz = 3*1.4;
cocosz = 7 *1.4;
%}

ustabc1 = [0.9290    0.6940    0.1250];
stabc1 = [0.4940    0.1840    0.5560];
ustabc2 = [0.4660    0.6740    0.1880];
stabc2 = [0.3010    0.7450    0.9330];
%% Prepare SSMdata
sSample = 1;
numoutdof = numel(outdof);

epsfssm = [];
for i = 1:numel(mus)
   mui = mus(i) * ones(1,numel(FRCom(i).FRC));
   %FRCom(i).FRC.Aout = FRCom(i).FRC.Aout(:,1)
   %FRCom(i).FRC.mus = mui;
   epsfssm = [epsfssm,mui];
end
FRCom = vertcat(FRCom.FRC);
stab  = [FRCom.stability];
Aout = [FRCom.Aout];
Aout = Aout(outdof:2:end);
omssm = [FRCom.Omega];


%% Prepare COCOdata
nSample = 1;
inSample = 4;
ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
    ampNames{k} = strcat('amp',num2str(outdof(k)));
end

omegacc = [];
epsfcc  = [];
stabcc = logical([]);
amp = cell(numoutdof,1);

for i=1:numel(bds)
    bd = bds(i).bd;
    bd = bd{1};
    omegai = coco_bd_col(bd, 'omega');
    epsfi = coco_bd_col(bd, 'eps');
    epsfi = ones(1,numel(epsfi))*mus(i);
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
    
    %% SSM Plot
    Aouti = Aout(i,:);
    
    oms = omssm(stab); eps = epsfssm(stab); am = Aouti(stab);
    plot3(oms(1:sSample:end),  eps(1:sSample:end),am(1:sSample:end),'bo','MarkerSize',ssmsz,'MarkerFaceColor','b');
    oms  = omssm(~stab); eps = epsfssm(~stab); am = Aouti(~stab);
    plot3(oms(1:sSample:end),  eps(1:sSample:end),am(1:sSample:end),'ro','MarkerSize',ssmsz,'MarkerFaceColor','r');
    
    if any(~stab)
    Leg{end+1} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable');
    end
    Leg{end+1} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable');
        
    %% Coco Plot
    %{
    om = omegacc(stabcc); am = amp{k}(stabcc); eps = epsfcc(stabcc);
    plot3(om(1:nSample:end),  eps(1:nSample:end),am(1:nSample:end), 'g*','MarkerSize', cocosz);%,...
        %'color',[0 0 0] );
    om = omegacc(~stabcc); am = amp{k}(~stabcc);eps = epsfcc(~stabcc);
    plot3(om(1:nSample:end), eps(1:nSample:end),am(1:nSample:end), 'c*', 'MarkerSize', cocosz);%, ...
        %'color',[0.9294    0.6941    0.1255] );
    %}
    % {   
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
        
    
    %% Labels and Legend
    legend(Leg{2:end},...
        'Interpreter','latex','Location','best');
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\mu$','interpreter','latex','FontSize',16);
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

