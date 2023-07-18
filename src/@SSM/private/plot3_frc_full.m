function plot3_frc_full(omega,epsilon,Znorm,outdof,Aout,stab,order,varargin)

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
    figure();view(3);hold on 
    plot3(omega(stab),epsilon(stab),Znorm(stab),'ob','MarkerSize',4) 
    plot3(omega(~stab),epsilon(~stab),Znorm(~stab),'or','MarkerSize',4)
     
    xlabel('$\Omega$','Interpreter','latex'); 
    ylabel('$\epsilon$','Interpreter','latex'); 
    zlabel('$\|\mathbf{z}\|_{L_2}$','Interpreter','latex'); 
    if all(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            'Interpreter','latex', 'Location','best','Box','off');
    elseif ~any(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
            'Interpreter','latex', 'Location','best','Box','off');    
    else
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
            'Interpreter','latex', 'Location','best','Box','off');
    end
    set(gca,'FontSize',14);
    grid on, axis tight;view(3);
    
    figure; hold on
    numOutdof = numel(outdof);
    for k=1:numOutdof
       subplot(numOutdof,1,k);view(3);hold on 
       plot3(omega(stab),epsilon(stab),Aout(stab,k),'ob','MarkerSize',4);      
       plot3(omega(~stab),epsilon(~stab),Aout(~stab,k),'or','MarkerSize',4);
       zk = strcat('$||z_{',num2str(outdof(k)),'}||_{\infty}$');
       zlabel(zk,'Interpreter','latex');
       set(gca,'FontSize',14);
       grid on, axis tight

    
    xlabel('$\Omega$','Interpreter','latex');
    ylabel('$\epsilon$','Interpreter','latex'); 
    if all(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            'Interpreter','latex', 'Location','best','Box','off');
    elseif ~any(stab)
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
            'Interpreter','latex', 'Location','best','Box','off');    
    else
        legend(strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - stable'),...
            strcat('SSM-$$\mathcal{O}(',num2str(order),')$$ - unstable'),...
            'Interpreter','latex', 'Location','best','Box','off');
    end
    end
         

end
