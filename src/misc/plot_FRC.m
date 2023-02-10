function plot_FRC(FRC,outdof,order,ParName,plotStyle,figs,color)
%% PLOT_FRC_FULL This function plots forced response curve.
%
% PLOT_FRC_FULL(FRC,OUTDOF,ORDER,PARNAME,PLOTSTYLE,FIGS,COLOR)
%
% FRC:       a struct array
% OUTDOF:    plot of FRC for these dofs
% ORDER:     expansion order of autonomous SSM
% PARNAME:   freq/amp FRC with varied forcing frequency or amplitude
% PLOTSTYLE: ploting with lines or circles
% COLOR:     color of lines/circles
switch ParName
    case 'freq'
        xlab = '$\Omega$';
        if isfield(FRC,'Omega')
            Par  = [FRC.Omega];
        elseif isfield(FRC,'om')
            Par  = [FRC.om];
        end
    case 'amp'
        xlab = '$\epsilon$';
        Par  = [FRC.epsilon];
    otherwise
        error('The ParName should be amp/freq');
end

% extract other relevant data
Aout  = [FRC.Aout];
stab  = [FRC.stability];
Znorm = [FRC.Znorm];
if isa(outdof,'function_handle');outdof = 1:round(numel(Aout)/numel(stab));end
numOutdof = numel(outdof);
numPts    = numel(stab);
Aout = reshape(Aout, [numOutdof, numPts]);
Aout = Aout';

switch plotStyle
    case 'lines'
        % plot by solid and dashed lines (solid/dashed: stable/unstable)
        % plot for znorm
        figure(figs(1)); hold on
        stab_plot(Par,Znorm,stab,order,color);
        
        % plot of special points (SN and HB)
        [SNHB,SNidx,HBidx] = get_SNHB(numPts,FRC);
        if SNHB
            plot_SNHB(Par(SNidx),Znorm(SNidx),Par(HBidx),Znorm(HBidx))
        end
        add_labels('$\Omega$','$\|\mathbf{z}\|_{L_2}$')
        
        % plot for Aout at outdofs
        figure(figs(2)); hold on
        for k=1:numOutdof
            if numOutdof > 1
                subplot(numOutdof,1,k); hold on
            end
            stab_plot(Par,Aout(:,k),stab,order,color);
            % plot of special points
            if SNHB
                plot_SNHB(Par(SNidx),Aout(SNidx,k),Par(HBidx),Aout(HBidx,k))
            end
            add_labels('$\Omega$',strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$'))
        end
        xlabel(xlab,'Interpreter','latex');
        
    case 'circles'
        % plot by blue and red cycles (blue/red: stable/unstable)
        % plot for znorm
        figure(figs(1)); hold on
        plot(Par(stab),Znorm(stab),'o','Color', color,'MarkerSize',10,'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'));
        plot(Par(~stab),Znorm(~stab),'s','Color', color, 'MarkerSize',10,'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'));
        add_labels('$\Omega$','$\|\mathbf{z}\|_{L_2}$')
        lgd = legend();
        set(lgd,'Interpreter','latex','Location','best');
        
        % plots for Aout at outdofs
        figure(figs(2)); hold on
        for k=1:numOutdof
            if numOutdof > 1
                subplot(numOutdof,1,k); hold on
            end
            plot(Par(stab),Aout(stab,k),'o','Color', color,'MarkerSize',10,'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'));
            plot(Par(~stab),Aout(~stab,k),'s','Color', color,'MarkerSize',10,'DisplayName',strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'));
            add_labels('$\Omega$',strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$'))
            lgd = legend();
            set(lgd,'Interpreter','latex','Location','best');
            xlabel('$\Omega$','Interpreter','latex');
        end
end

end


function add_legends(stab,order)
Leg0 =  get(legend(),'string');
if all(stab)
    legend(Leg0{:},strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
        'Interpreter','latex','Location','best');
elseif ~any(stab)
    legend(Leg0{:},strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
        'Interpreter','latex','Location','best');
else
    legend(Leg0{:},strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
        strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
        'Interpreter','latex','Location','best');
end

end

function [SNHB,SNidx,HBidx] = get_SNHB(numPts,FRC)
% get indices of bifurcation points on the FRC
SNidx = false(numPts,1);
HBidx = false(numPts,1);
if isfield(FRC,'isSN')
    SNidx = [FRC.isSN];
end
if isfield(FRC,'isHB')
    HBidx = [FRC.isHB];
end
SNHB = any([SNidx(:);HBidx(:)]);
end

function plot_SNHB(xSN,ySN,xHB,yHB)
thm = struct();
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
    'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
    'black', 'MarkerFaceColor', 'white'};
SNfig = plot(xSN, ySN,thm.SN{:});
set(get(get(SNfig,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
HBfig = plot(xHB,yHB,thm.HB{:});
set(get(get(HBfig,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');
end

function add_labels(xlab,ylab)
xlabel(xlab,'Interpreter','latex');
ylabel(ylab,'Interpreter','latex');
set(gca,'FontSize',14);
grid on, axis tight; legend boxoff;
end