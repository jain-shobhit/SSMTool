function plot3_frc_full(omega,epsilon,Znorm,outdof,Aout,stab,order,varargin)

if isa(outdof,'function_handle');outdof = 1:size(Aout,2);end
if numel(varargin)>0 && strcmp(varargin{1},'lines')
    % plot by solid and dashed lines (solid/dashed: stable/unstable)
    if numel(varargin)>1 && strcmp(varargin{2},'BP')
        figure(gcf); hold on
        numOutdof = numel(outdof);
        if numOutdof>1
            for k=1:numOutdof
               subplot(numOutdof,1,k); hold on
               stab_plot3(omega,epsilon,Aout(:,k),stab,order,'nolegends');
            end
        else
            stab_plot3(omega,epsilon,Aout,stab,order,'nolegends');
        end
    else
        figure; hold on
        if isempty(stab)
            plot3(omega,epsilon,Znorm,varargin{2},'LineWidth',2);
            leg = strcat('SSM-$$\mathcal{O}(',num2str(order),')$$');
            legend(leg,'Interpreter','latex');
        else
            stab_plot3(omega,epsilon,Znorm,stab,order);
        end
        xlabel('$\Omega$','Interpreter','latex'); 
        ylabel('$\epsilon$','Interpreter','latex');
        zlabel('$\|\mathbf{z}\|_{L_2}$','Interpreter','latex'); 
        set(gca,'FontSize',14);
        view([50 15]);
        grid on, axis tight; legend boxoff; box on
        figure; hold on
        numOutdof = numel(outdof);
        for k=1:numOutdof
           subplot(numOutdof,1,k); hold on
           if isempty(stab)
               plot3(omega,epsilon,Aout(:,k),varargin{2},'LineWidth',2);
               legend(leg,'Interpreter','latex');
           else
               stab_plot3(omega,epsilon,Aout(:,k),stab,order);
           end
           zk = strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$');
           zlabel(zk,'Interpreter','latex');
           set(gca,'FontSize',14);
           view([50 15]);
           grid on, axis tight; legend boxoff; box on
           xlabel('$\Omega$','Interpreter','latex'); 
           ylabel('$\epsilon$','Interpreter','latex');
        end    
    end
else
    % plot by blue and red cycles (blue/red: stable/unstable)
    % BUG <NEED TO BE FIXED FOR THE CASE OF BP PLOTS>
    figure; hold on
    plot3(omega(stab),epsilon(stab),Znorm(stab),'ob','MarkerSize',10)
    plot3(omega(~stab),epsilon(stab),Znorm(~stab),'or','MarkerSize',10)
    xlabel('$\Omega$','Interpreter','latex'); 
    ylabel('$\epsilon$','Interpreter','latex'); 
    zlabel('$\|\mathbf{z}\|_{L_2}$','Interpreter','latex'); 
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
    grid on, axis tight;
    figure; hold on
    numOutdof = numel(outdof);
    for k=1:numOutdof
       subplot(numOutdof,1,k); hold on
       plot3(omega(stab),epsilon(stab),Aout(stab,k),'ob','MarkerSize',10);
       plot3(omega(~stab),epsilon(~stab),Aout(~stab,k),'or','MarkerSize',10);
       zk = strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$');
       zlabel(zk,'Interpreter','latex');
       set(gca,'FontSize',14);
       grid on, axis tight
    end
    xlabel('$\Omega$','Interpreter','latex'); 
    ylabel('$\epsilon$','Interpreter','latex'); 
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
end
end