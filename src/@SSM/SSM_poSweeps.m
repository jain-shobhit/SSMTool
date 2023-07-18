function FRC = SSM_poSweeps(obj,oid,modes,order,epSamp,omRange,varargin)
% SSM_POSWEEPS This function performs a family of continuations of 
% periodic orbits of reduced dynamics. Continuation is performed with varied
% forcing frequency for each forcing amplitude.
% The continuation here starts from the guess of an initial solution. 
%
% FRC = SSM_POSWEEPS(OBJ,OID,MODES,ORDER,EPSAMPS,OMRANGE, VARARGIN)
%
% oid    :      runid of continuation
% modes  :      master spectral subspace
% order  :      expansion order of SSM
% epSamp :      sampled forcing amplitudes for forced response curves
% omRange:      continuation domain of forcing frequency, which should be near the
%               value of natural frequency with index 1
% varargin:     parameter that can be used to suppres the plot of FRCs
%
% FRC:          Forced Response Curve data struct
%
% See also: COCO_POSWEEPS, SSM_LVLSWEEPS, SSM_EPSWEEPS

if~isempty( varargin)
    plotfrc = varargin{1}; 
end

outdof = obj.FRCOptions.outdof;
FRCom  = cell(numel(epSamp),1);

ii = 1;
for ep = epSamp
    oideps = strcat(oid,num2str(ep));
    if obj.System.order == 1
        obj.System.Fext.epsilon = ep;
    else
        obj.System.fext.epsilon = ep;       
    end

    FRCom{ii} = FRC_cont_po(obj,oideps,modes,order,omRange);
    ii = ii +1;
end

%% results extraction

FRCom = horzcat(FRCom{:});

% save results
FRC = struct();
FRC.SSMorder   = order;
FRC.SSMnonAuto = obj.Options.contribNonAuto;
FRC.FRCom = FRCom;
FRC.eps = epSamp;



%% Visualization of results
if~isempty( varargin)
    if plotfrc
        sweeps_plot(numel(FRCom),FRCom,outdof,order);
    end
else
   sweeps_plot(numel(FRCom),FRCom,outdof,order);
end
end

function sweeps_plot(numLabs,FRCom,outdof,order)
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
    % 2D plot
    figure; hold on
    for k=1:numLabs
        om = FRCom{k}.om;
        Aout = FRCom{k}.Aout(i,:);
        stab = FRCom{k}.stability;
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
        %HBidx = FRCom{k}.HBidx;
        %HBfig = plot(om(HBidx),Aout(HBidx),thm.HB{:});
        %set(get(get(HBfig,'Annotation'),'LegendInformation'),...
        %    'IconDisplayStyle','off');         
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
        Aout = FRCom{k}.Aout(i,:);
        stab = FRCom{k}.stability;
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
        %HBidx = FRCom{k}.HBidx;
        %HBfig = plot3(om(HBidx),epsf(HBidx),Aout(HBidx),thm.HB{:});
        %set(get(get(HBfig,'Annotation'),'LegendInformation'),...
        %    'IconDisplayStyle','off');        
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

function stab_plot3(x,y,z,S,order,varargin)
% This function is adapted from stab_plot in coco_plot_bd in coco

II  = 1;
EI  = numel(x);
I   = II;
figs = [];
stab = [];
ST = cell(2,1);
%ST{1} = {'ro','MarkerSize', 4};
%ST{2} = {'bo','MarkerSize', 4};
ST{1} = {'b--','LineWidth',1.5};
ST{2} = {'b-','LineWidth',1.5};
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
    Leg{1} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable');
    Leg{2} = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable');
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

function [I, I0] = next_index(EI, S, I, II)
while (I<EI && S(II)==S(I)); I=I+1; end
I0 = I;
end
