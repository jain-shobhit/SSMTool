function bds = coco_poSweeps(obj,epSamp,omRange,varargin)
% COCO_POSWEEPS This function performs continuation of a family of periodic
% orbits with varying frequency for a given set of sampling forcing
% amplitudes. For each amplitude the FRC is computed.
%
% BDS = COCO_POSWEEPS(OBJ,EPSAMP,OMRANGE, varargin)
%
% obj:      instance of DS class
% epSamp:   Sampling forcing amplitudes
% omRange:  continuation domain of forcing frequency
% varargin: parameter that can be used to suppress the plotting of results
%
% bds:      bifurcation data struct that contains the information on the
%           continuation problems that are investigated
%
% See also: SSM_POSWEEPS

if~isempty( varargin)
    plotfrc = varargin{1}; 
end
bds = cell(1,numel(epSamp));
ii = 1;
for eps = epSamp  
    
    % Set current epsilon
    if obj.system.order == 1
        obj.system.Fext.epsilon = eps;
    else
        obj.system.fext.epsilon = eps;
    end
    
    % Reset the initial forcing frequency
    obj.system.Omega = omRange(1);
    
    % Get individual FRCs
    bds{ii} = obj.extract_FRC( omRange,false); %dont plot frcs 
    ii = ii +1;
end

% Prepare output
bds = vertcat(bds{:});

if~isempty( varargin)
    if plotfrc
        Sweepplot( bds,obj.outdof);
    end
else
    Sweepplot( bds,obj.outdof);
    
end

end

function Sweepplot( bds,outdof)
%% Plot Setting
cocosz = 5 *1.4;

%% Prepare COCOdata
numoutdof = numel(outdof);
nSample = 1;

ampNames = cell(1, numel(numoutdof));
for k = 1:numel(outdof)
    ampNames{k} = strcat('amp',num2str(outdof(k)));
end

omegacc = [];
epsfcc  = [];
stabcc = logical([]);
amp = cell(numoutdof,1);


for k=1:numel(bds)
    bd = bds{k};
    
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

figure(gcf)
hold on
if numoutdof > 1
    for i=1:numoutdof
        subplot(numoutdof,1,i); hold on
        
        %% Coco Plot
        
        om = omegacc(~stabcc); am = amp{i}(~stabcc);eps = epsfcc(~stabcc);
        plot3(om(1:nSample:end), eps(1:nSample:end),am(1:nSample:end), 'c*', 'MarkerSize', cocosz,'DisplayName','collocation - unstable')%, 'color',[0.9294    0.6941    0.1255] );
        om = omegacc(stabcc); am = amp{i}(stabcc); eps = epsfcc(stabcc);
        plot3(om(1:nSample:end),  eps(1:nSample:end),am(1:nSample:end), 'g*','MarkerSize', cocosz,'DisplayName','collocation - stable')%,'color',[0 0 0] );
    end
else
    om = omegacc(~stabcc); am = amp{1}(~stabcc);eps = epsfcc(~stabcc);
    plot3(om(1:nSample:end), eps(1:nSample:end),am(1:nSample:end), 'c*', 'MarkerSize', cocosz,'DisplayName','collocation - unstable')%, 'color',[0.9294    0.6941    0.1255] );
    om = omegacc(stabcc); am = amp{1}(stabcc); eps = epsfcc(stabcc);
    plot3(om(1:nSample:end),  eps(1:nSample:end),am(1:nSample:end), 'g*','MarkerSize', cocosz,'DisplayName','collocation - stable')%,'color',[0 0 0] );
    
end
    %% Labels and Legend
    % legend(Leg{2:end},  'Interpreter','latex','Location','best');
    % xlabel('$\Omega$','interpreter','latex','FontSize',16);
    % ylabel('$\mu$','interpreter','latex','FontSize',16);
    % zk = strcat('$||z_{',num2str(outdof(i)),'}||_{\infty}$');
    % zlabel(zk,'Interpreter','latex');
    %% Parameters for figure size, label size, legend size
    % set(gca,'fontsize',axis_size)
    
    %% Size of window
    % set(gcf, 'Position', [0 0 xwidth ywidth])
    
    % grid on, axis tight;legend boxoff;
    % view([-1,-1,1]);

end

