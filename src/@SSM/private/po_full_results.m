function FRC = po_full_results(runid,sampStyle,ispolar,ampNames,saveIC,branch,varargin)
% PO_FULL_RESULTS This function extract results of reduced dynamics at
% sampled forcing frequencies/amplitudes in continuation of periodic orbits
% The response is computed in the full system
%
% FULLSAMP =
% PO_REDUCED_RESULTS(RUNID,SAMPSTYLE,ISPOLAR,ISOMEGA,AMPNAMES,SAVEIC,BRANCH,VARARGIN)
%
% runid:     id of continuation run of equilibrium points
% sampStyle: samping style of excitation frequencies
% ispolar:   coordinate representation of reduced dynamics
% ampNames:  names of the amplitudes that are output to FRC computation
% saveIC:    whether to save initial condition of FRC orbits in full
%            coordinates
% branch:    if branch = primary   -> primary branch
%            if branch = secondary -> append results from secondary branch
% varargin:  FRC struct array of primary branch

%% extract results of reduced dynamics at sampled frequency/forcing
bd = coco_bd_read(runid);
m  = numel(ampNames)-1; % last one is Znorm 

Aout_frc = [];
Znorm_frc = [];

if strcmp(sampStyle,'cocoBD')
    om  = coco_bd_col(bd, 'om');
    
    if ispolar
        % Not implemented yet
    else
        for k=1:m    
            % Read amplitudes of outdofs
            Aout_k =  coco_bd_col(bd, ampNames{k});            
            Aout_frc = [Aout_frc; Aout_k];
           
        end
        Znorm_frc = coco_bd_col(bd, 'Znorm'); 
        if saveIC
        Z0_frc = coco_bd_col(bd, 'IC');
        end
    end
    epsf = coco_bd_col(bd, 'eps');
    

    % Floquet multiplier give stability
    stab  = ~coco_bd_col(bd, 'po.test.USTAB');

    SNidx = coco_bd_idxs(bd, 'SN');
    PDidx = coco_bd_idxs(bd, 'PD');
    TRidx = coco_bd_idxs(bd, 'TR');
    BPidx = coco_bd_idxs(bd, 'BP');

else
    if strcmp(sampStyle, 'uniform')
        error('Sampling Style "uniform" not supported for po');
    elseif strcmp(omegaSampStyle, 'cocoOut')
        error('Sampling Style "cocoOut" not supported for po');
    end

end
%% Append if results come from secondary branch

switch branch
    case 'primary'
        FRC = struct();
        FRC.om  = om;
        FRC.Aout  = Aout_frc;
        FRC.Zout  = []; % Not stored for po
        FRC.Znorm = Znorm_frc;
        FRC.ep  = epsf;
        FRC.stability  = stab;
        FRC.SNidx = SNidx;
        FRC.PDidx = PDidx;
        FRC.TRidx = TRidx;
        FRC.BPidx = BPidx;
        
        if saveIC
            FRC.Z0_frc    = Z0_frc; % initial state
        else
            FRC.Z0_frc    = []; % initial state
            
        end
        
    case 'secondary'
        FRC = varargin{1};
        FRC.om        = [FRC.om;om];
        FRC.Aout  = [FRC.Aout_frc;Aout_frc];
        FRC.Zout  = []; % Not stored for po
        FRC.Znorm = [FRC.Znorm_frc;Znorm_frc];
        FRC.ep        = [FRC.ep;epsf];
        FRC.stability    = [FRC.st;stab];
        FRC.SNidx = [FRC.SNidx,SNidx];
        FRC.PDidx = [FRC.PDidx,PDidx];
        FRC.TRidx = [FRC.TRidx,TRidx];
        FRC.BPidx = [FRC.BPidx,BPidx];
        
        if saveIC
            FRC.Z0_frc    = [FRC.Z0_frc;Z0_frc]; % initial state
        end
end
%{
if numel(varargin)==1
    
%% plot continuation results in normal coordinates
thm = struct( 'special', {{'SN', 'HB'}});
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.HB = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
figure;
for k=1:m
    subplot(m,1,k)
    if ispolar
        rhok = strcat('rho',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', rhok, @(x) abs(x)); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', rhok, @(x) abs(x)); hold on
            else
                coco_plot_bd(thm, runid, 'eps', rhok, @(x) abs(x)); hold on
            end
        end
    else
        Rezk = strcat('Rez',num2str(k));
        Imzk = strcat('Imz',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', {Rezk,Imzk}, @(x,y) sqrt(x.^2+y.^2)); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', {Rezk,Imzk}, @(x,y) sqrt(x.^2+y.^2)); hold on
            else
                coco_plot_bd(thm, runid, 'eps', {Rezk,Imzk}, @(x,y) sqrt(x.^2+y.^2)); hold on
            end
        end
    end
    grid on; box on; 
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    if isempty(isomega)
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$\epsilon$','interpreter','latex','FontSize',16);
        zlabel(strcat('$\rho_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    else
        if isomega
            xlabel('$\Omega$','interpreter','latex','FontSize',16);
        else
            xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        end
        ylabel(strcat('$\rho_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    end
end


figure;
for k=1:m
    subplot(m,1,k)
    if ispolar
        rhok = strcat('rho',num2str(k));
        thk  = strcat('th',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', {thk, rhok}, @(x,y) x+0.5*(1-sign(y))*pi); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', {thk, rhok}, @(x,y) x+0.5*(1-sign(y))*pi); hold on
            else
                coco_plot_bd(thm, runid, 'eps', {thk, rhok}, @(x,y) x+0.5*(1-sign(y))*pi); hold on
            end
        end
    else
        Rezk = strcat('Rez',num2str(k));
        Imzk = strcat('Imz',num2str(k));
        if isempty(isomega)
            coco_plot_bd(thm, runid, 'om', 'eps', {Rezk,Imzk}, @(x,y) atan2(y,x)); hold on
        else
            if isomega
                coco_plot_bd(thm, runid, 'om', {Rezk,Imzk}, @(x,y) atan2(y,x)); hold on
            else
                coco_plot_bd(thm, runid, 'eps', {Rezk,Imzk}, @(x,y) atan2(y,x)); hold on
            end
        end
    end
    grid on; box on; 
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    if isempty(isomega)
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$\epsilon$','interpreter','latex','FontSize',16);
        zlabel(strcat('$\theta_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    else
        if isomega
            xlabel('$\Omega$','interpreter','latex','FontSize',16);
        else
            xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        end
        ylabel(strcat('$\theta_',num2str(k),'$'),'interpreter','latex','FontSize',16);
    end
end

end
%}
end
