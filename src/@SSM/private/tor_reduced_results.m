function redSamp = tor_reduced_results(runid,ispolar,isomega,mFreqs,nt)
% PO_REDUCED_RESULTS This function extract results of reduced dynamics at
% sampled forcing frequencies/amplitudes in continuation of periodic orbits
%
% runid:     id of continuation run of equilibrium points
% ispolar:   coordinate representation of reduced dynamics
% isomega:   is the continuation parameter forcing frequency
% mFreqs:    internal resonance relation vector
% nt:        number of time points for discretizing trajectory

%% extract results of reduced dynamics at sampled frequency/forcing
bd   = coco_bd_read(runid);
labs = coco_bd_labs(bd,'all');
nlab = numel(labs);
sols = cell(nlab,1);
omeg = zeros(nlab,1);
epsf = zeros(nlab,1);
for k = 1:nlab
    sol = tor_read_solution('',runid, labs(k));
    omeg(k) = sol.p(1);
    epsf(k) = sol.p(2);
    sols{k} = sol;
end

redSamp = struct();
redSamp.om  = omeg;
redSamp.ep  = epsf;
redSamp.lab = labs;

%% plot continuation path in normal coordinates
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
        coco_plot_bd(runid, 'om', 'om1');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$\omega_1$','interpreter','latex','FontSize',16);
        subplot(2,1,2);
        coco_plot_bd(runid, 'om', 'om2');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\Omega$','interpreter','latex','FontSize',16);
        ylabel('$\omega_2$','interpreter','latex','FontSize',16);           
    else
        figure; 
        subplot(2,1,1);
        coco_plot_bd(runid, 'eps', 'om1');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        ylabel('$\omega_1$','interpreter','latex','FontSize',16);
        subplot(2,1,2);
        coco_plot_bd(runid, 'eps', 'om2');
        grid on; box on; 
        set(gca,'LineWidth',1.2);
        set(gca,'FontSize',12);
        xlabel('$\epsilon$','interpreter','latex','FontSize',16);
        ylabel('$\omega_2$','interpreter','latex','FontSize',16);
    end 
end

%% construct torus in reduced system using interpolation of periodic cycle
disp('Constructing 3D torus in reduced dynamical system');
qTr  = cell(nlab,1);
tTr  = cell(nlab,1);
nSeg = zeros(nlab,1);
angs = linspace(0,2*pi,100);
for i=1:numel(labs)
    soli = sols{i};
    tbp = soli.tbp;
    xbp = soli.xbp; % xbp = zeros(nt, data.bc_data.dim, N+1);
    % recover to high fidelty based on Fourier expansion
    xbp = tor_interpolation(xbp,angs);
    om  = soli.p(1);
    assert(abs(om-omeg(i))<1e-3*omeg(i), 'Read wrong solution from data');
    fprintf('Interpolation at frequency %d\n', om);
    ratio = ceil(soli.p(1)/soli.p(end-1));
    tsamp = linspace(0,tbp(end),ratio*nt);
    [~,dim, numSegs] = size(xbp);
    ys = zeros(ratio*nt,dim,numSegs);
    for k=1:numSegs
        % Interpolation at fine mesh
        xk  = interp1(tbp, xbp(:,:,k), tsamp, 'pchip'); % Update basepoint values
        % Real to complex
        if ispolar
            zk = xk(:,1:2:end-1).*exp(1j*xk(:,2:2:end));
        else
            zk = xk(:,1:2:end-1)+1j*xk(:,2:2:end);
        end
        % \mathbb{S}^2 --> \mathbb{S}^3, slow-time --> regular time
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

end
