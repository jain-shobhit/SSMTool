function plot_frc_full(Par,Znorm,outdof,Aout,stab,order,ParName,varargin)

switch ParName
    case 'freq'
        xlab = '$\Omega$';
    case 'amp'
        xlab = '$\epsilon$';
    otherwise
        error('The ParName should be amp/freq');
end

if isa(outdof,'function_handle');outdof = 1:size(Aout,2);end

if numel(varargin)>0 && strcmp(varargin{1},'lines')
    % plot by solid and dashed lines (solid/dashed: stable/unstable)    
    if numel(varargin)>1 && ischar(varargin{2})
        if strcmp(varargin{2},'BP')
        figure(gcf); hold on
%         stab_plot(Par,Znorm,stab,order,'nolegends');
        numOutdof = numel(outdof);
        if numOutdof>1
            for k=1:numOutdof
               subplot(numOutdof,1,k); hold on
               stab_plot(Par,Aout(:,k),stab,order,'blue','nolegends');
            end
        else
            stab_plot(Par,Aout,stab,order,'blue','nolegends');
        end
        end
    else
        figure; hold on
        stab_plot(Par,Znorm,stab,order,'blue');
        SNHB = false;
        if numel(varargin)>1 && iscell(varargin{2})
            SNHB = true;
            thm = struct();
            thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
              'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
              'cyan', 'MarkerFaceColor', 'white'};
            thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
              'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
              'black', 'MarkerFaceColor', 'white'};
            SNidx = varargin{2}{1};
            HBidx = varargin{2}{2};
            SNfig = plot(Par(SNidx),Znorm(SNidx),thm.SN{:});
            set(get(get(SNfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
            HBfig = plot(Par(HBidx),Znorm(HBidx),thm.HB{:});
            set(get(get(HBfig,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');            
        end
        xlabel(xlab,'Interpreter','latex'); 
        ylabel('$\|\mathbf{z}\|_{L_2}$','Interpreter','latex'); 
        set(gca,'FontSize',14);
        grid on, axis tight; legend boxoff;
        figure; hold on
        numOutdof = numel(outdof);
        for k=1:numOutdof
           subplot(numOutdof,1,k); hold on
           stab_plot(Par,Aout(:,k),stab,order,'blue');
           if SNHB
                SNfig = plot(Par(SNidx),Aout(SNidx,k),thm.SN{:});
                set(get(get(SNfig,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');
                HBfig = plot(Par(HBidx),Aout(HBidx,k),thm.HB{:});
                set(get(get(HBfig,'Annotation'),'LegendInformation'),...
                'IconDisplayStyle','off');            
           end               
           zk = strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$');
           ylabel(zk,'Interpreter','latex');
           set(gca,'FontSize',14);
           grid on, axis tight; legend boxoff;
        end
        xlabel(xlab,'Interpreter','latex');    
    end
else
    % plot by blue and red cycles (blue/red: stable/unstable)
    % BUG <NEED TO BE FIXED FOR THE CASE OF BP PLOTS - not used though>
    figure; hold on
    plot(Par(stab),Znorm(stab),'ob','MarkerSize',10)
    plot(Par(~stab),Znorm(~stab),'or','MarkerSize',10)
    xlabel('$\Omega$','Interpreter','latex'); 
    ylabel('$\|\mathbf{z}\|_{L_2}$','Interpreter','latex'); 
    if all(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            'Interpreter','latex');
    elseif ~any(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
            'Interpreter','latex');    
    else
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),'Interpreter','latex');
    end
    set(gca,'FontSize',14);
    grid on, axis tight; legend boxoff;
    figure; hold on
    numOutdof = numel(outdof);
    for k=1:numOutdof
       subplot(numOutdof,1,k); hold on
       plot(Par(stab),Aout(stab,k),'ob','MarkerSize',10);
       plot(Par(~stab),Aout(~stab,k),'or','MarkerSize',10);
       zk = strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$');
       ylabel(zk,'Interpreter','latex');
       set(gca,'FontSize',14);
       grid on, axis tight
    end
    xlabel('$\Omega$','Interpreter','latex'); 
    if all(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            'Interpreter','latex');
    elseif ~any(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable',toolbox,')'),...
            'Interpreter','latex');    
    else
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),'Interpreter','latex');
    end
    legend boxoff;
end
end
