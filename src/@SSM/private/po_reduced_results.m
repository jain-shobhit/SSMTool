function redSamp = po_reduced_results(runid,ispolar,isomega,mFreqs,nt,isuniform)
% PO_REDUCED_RESULTS This function extract results of reduced dynamics at
% sampled forcing frequencies/amplitudes in continuation of periodic orbits
%
% runid:     id of continuation run of equilibrium points
% ispolar:   coordinate representation of reduced dynamics
% isomega:   is the continuation parameter forcing frequency
% mFreqs:    internal resonance relation vector
% nt:        number of time points for discretizing trajectory
% isuniform: is sampling style of omega/epsilon 'uniform'

%% extract results of reduced dynamics at sampled frequency/forcing
bd   = coco_bd_read(runid);
if isuniform
    labs = coco_bd_labs(bd, 'UZ');
else
    labs = coco_bd_labs(bd,'all');
end
nlab = numel(labs);
sols = cell(nlab,1);
if isempty(isomega)
    stab = nan;
else
    stab = false(nlab,1);
end
omeg = zeros(nlab,1);
epsf = zeros(nlab,1);
for k = 1:nlab
    sol = po_read_solution('',runid, labs(k));
    if ~isempty(isomega)
        stab(k) = all(abs(sol.po_test.la)<1);
    end
    omeg(k) = sol.p(1);
    epsf(k) = sol.p(2);
    sols{k} = sol;
end

redSamp = struct();
redSamp.om  = omeg;
redSamp.st  = stab;
redSamp.ep  = epsf;
redSamp.lab = labs;

%% plot continuation path in normal coordinates
thm = struct( 'special', {{'SN', 'TR', 'PD'}});
thm.SN = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'cyan', 'Marker', 'o', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'cyan', 'MarkerFaceColor', 'white'};
thm.TR = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 's', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
thm.PD = {'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'red', 'Marker', 'd', 'MarkerSize', 8, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'white'};
if isempty(isomega)
    figure;
    coco_plot_bd(runid, 'om', 'eps');
    grid on; box on; 
    set(gca,'LineWidth',1.2);
    set(gca,'FontSize',14);
    xlabel('$\Omega$','interpreter','latex','FontSize',16);
    ylabel('$\epsilon$','interpreter','latex','FontSize',16);
else
    if isomega
        figure; 
        subplot(2,1,1);
        coco_plot_bd(thm, runid, 'om', 'po.period');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$T$','interpreter','latex','FontSize',16);
        subplot(2,1,2);
        coco_plot_bd(thm, runid, 'om', '||x||_{2,D}');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$||x-\bar{x}||_{\mathcal{L}_2[0,1]}$','interpreter','latex','FontSize',16);           
    else
        figure; 
        subplot(2,1,1);
        coco_plot_bd(thm, runid, 'eps', 'po.period');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        ylabel('$T$','interpreter','latex','FontSize',16);
        subplot(2,1,2);
        coco_plot_bd(thm, runid, 'eps', '||x||_{2,D}');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        ylabel('$||x-\bar{x}||_{\mathcal{L}_2[0,1]}$','interpreter','latex','FontSize',16);
    end 
end

%% construct torus in reduced system using interpolation of periodic cycle
disp('Constructing torus in reduced dynamical system');
qTr  = cell(nlab,1);
tTr  = cell(nlab,1);
nSeg = zeros(nlab,1);
dim  = size(sol.xbp,2);
for i=1:numel(labs)
    soli = sols{i};
    tbp = soli.tbp;
    xbp = soli.xbp;
    om  = soli.p(1);
    Ti  = soli.T;
    assert(abs(om-omeg(i))<1e-3*omeg(i), 'Read wrong solution from data');
    fprintf('Interpolation at (omega,epsilon) = (%d,%d)\n', om, soli.p(2));
    tsamp   = linspace(0,2*pi/om/min(mFreqs),nt);
    numSegs = numel(tbp);
    ys = zeros(nt,dim,numSegs);
    for k=1:numSegs
        % shift tsamp such that the initial time is the same as
        % tbp(k)
        tsampk = tsamp+tbp(k);
        tsampk = mod(tsampk,Ti); % mod with the period of the periodic cycle        
        xk  = interp1(tbp, xbp, tsampk, 'pchip'); % Update basepoint values
        if ispolar
            zk = xk(:,1:2:end-1).*exp(1j*xk(:,2:2:end));
        else
            zk = xk(:,1:2:end-1)+1j*xk(:,2:2:end);
        end
        zk  = zk.*exp(1j*(om*mFreqs.*tsamp'));
        ys(:,1:2:end-1,k) = real(zk);
        ys(:,2:2:end,k)   = imag(zk);
    end
    qTr{i} = ys;
    tTr{i} = tsamp;
    nSeg(i) = numSegs;               
end      
redSamp.qTr = qTr;
redSamp.tTr = tTr;
redSamp.nSeg = nSeg;

%% plot a sample of torus 
% The torus corresponds to the priodic orbit with maximmal deviation of 
% the time-rescaled periodic orbit from its state-space mean
disp('Illustration of the construction of torus in reduced dynamical system');
if isuniform
    idxpo  = coco_bd_idxs(bd, 'UZ');
else
    idxpo  = coco_bd_idxs(bd, 'all');
end
radius = coco_bd_col(bd, '||x||_{2,D}');
[~,idxMaxRadius]  = max(radius(idxpo));
labMaxRadius = labs(idxMaxRadius);
labMaxRadius = find(labMaxRadius==labs);
ys      = qTr{labMaxRadius};
numSegs = nSeg(labMaxRadius);
fprintf('Visualization of torus at (omega,epsilon)=(%d,%d)\n',...
    omeg(labMaxRadius),epsf(labMaxRadius));

% visualization of quasi-periodic trajectories
ya    = ys(1,:,:);
ya    = reshape(ya,[dim,numSegs]);
ytube = permute(ys,[3 1 2]);
figure; hold on
dnum = round(numSegs/64);
if dim<3
    % plot of x1-t-x2
    ts = tTr{labMaxRadius};
    ts = repmat(ts, [numSegs,1]); 
    h = surf(ytube(:,:,1),ts, ytube(:,:,2));
else
    h = surf(ytube(:,:,1),ytube(:,:,2),ytube(:,:,3));
end
set(h,'edgecolor','none')
% set(h, 'facecolor', [0.5 0.8 0.8])
set(h,'FaceAlpha',0.5);
view([1,1,1])
grid on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
if dim<3
    xlabel('$\mathrm{Re}(z_1)$','interpreter','latex','FontSize',16);
    ylabel('$t$','interpreter','latex','FontSize',16);
    zlabel('$\mathrm{Im}(z_1)$','interpreter','latex','FontSize',16);
    plot3(ya(1,:),ts(:,1)',ya(2,:),'b-','LineWidth',2);
    for k=1:dnum:numSegs 
        yk = ys(:,:,k);
        plot3(yk(1,1),ts(k,1),yk(1,2),'ko','LineWidth',3);
        plot3(yk(:,1),ts(k,:)',yk(:,2),'k-'); 
        plot3(yk(end,1),ts(k,end),yk(end,2),'bo','LineWidth',3);
        pause(0.2) 
    end    
else
    xlabel('$\mathrm{Re}(z_1)$','interpreter','latex','FontSize',16);
    ylabel('$\mathrm{Im}(z_1)$','interpreter','latex','FontSize',16);
    zlabel('$\mathrm{Re}(z_2)$','interpreter','latex','FontSize',16);
    plot3(ya(1,:),ya(2,:),ya(3,:),'b-','LineWidth',2);
    for k=1:dnum:numSegs 
        yk = ys(:,:,k);
        plot3(yk(1,1),yk(1,2),yk(1,3),'ko','LineWidth',3);
        plot3(yk(:,1),yk(:,2),yk(:,3),'k-'); 
        plot3(yk(end,1),yk(end,2),yk(end,3),'bo','LineWidth',3);
        pause(0.2) 
    end
end

% visualization of evoluation of closed curves
figure; hold on
if dim<3
    h = surf(ytube(:,:,1),ts, ytube(:,:,2));
else
    h = surf(ytube(:,:,1),ytube(:,:,2),ytube(:,:,3));
end
set(h,'edgecolor','none')
set(h,'FaceAlpha',0.5);
view([1,1,1])
grid on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
if dim<3
    xlabel('$\mathrm{Re}(z_1)$','interpreter','latex','FontSize',16);
    ylabel('$t$','interpreter','latex','FontSize',16);
    zlabel('$\mathrm{Im}(z_1)$','interpreter','latex','FontSize',16);

    for k=[1,round(nt/64):round(nt/64):nt-1, nt]
        yk = ys(k,:,:);
        yk = reshape(yk,[dim,numSegs]);
        if k==1 || k==nt
            plot3(yk(1,:),ts(:,k)',yk(2,:),'b-','LineWidth',2); pause(0.1)
        else
            plot3(yk(1,:),ts(:,k)',yk(2,:),'r-'); pause(0.1)
        end
    end
else
    xlabel('$\mathrm{Re}(z_1)$','interpreter','latex','FontSize',16);
    ylabel('$\mathrm{Im}(z_1)$','interpreter','latex','FontSize',16);
    zlabel('$\mathrm{Re}(z_2)$','interpreter','latex','FontSize',16);

    for k=[1,round(nt/64):round(nt/64):nt-1, nt]
        yk = ys(k,:,:);
        yk = reshape(yk,[dim,numSegs]);
        if k==1 || k==nt
            plot3(yk(1,:),yk(2,:),yk(3,:),'b-','LineWidth',2); pause(0.1)
        else
            plot3(yk(1,:),yk(2,:),yk(3,:),'r-'); pause(0.1)
        end
    end
end

end
