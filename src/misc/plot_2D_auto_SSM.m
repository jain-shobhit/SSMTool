function plot_2D_auto_SSM(W,rhosamp,plotdofs,varargin)
% PLOT_2D_AUTO_SSM This function returns the plot of 2 dimensional
% autonomous SSM. We first generate the grids of parameterization
% coordinates based on rhosamp X thetasamp(0:2*pi/128:2*pi). Then we map
% these grids to the full system based on the expansion W of SSM at
% plotdofs. Note that plotdofs should have 3 dofs.

assert(numel(plotdofs)==3,'the number of dofs is not three');
% generate grids
[RHO,THETA] = meshgrid(rhosamp,0:2*pi/128:2*pi);
zdof1 = zeros(size(RHO));
zdof2 = zeros(size(RHO));
zdof3 = zeros(size(RHO));
for k=1:129
    pk = RHO(k,:).*exp(1i*THETA(k,:));
    zk = reduced_to_full_traj([],[pk;conj(pk)],W);
    zk = zk(plotdofs,:);
    zdof1(k,:) = zk(1,:);
    zdof2(k,:) = zk(2,:);
    zdof3(k,:) = zk(3,:);
end

figure; hold on
h = surf(zdof1,zdof2,zdof3,'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 1.0, ...
    'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
    'LineWidth', 0.5);
view([1,1,1])
grid on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
if isempty(varargin)
    xlabel(['$z_{\mathrm{',num2str(plotdofs(1)),'}}$'],'interpreter','latex','FontSize',16);
    ylabel(['$z_{\mathrm{',num2str(plotdofs(2)),'}}$'],'interpreter','latex','FontSize',16);
    zlabel(['$z_{\mathrm{',num2str(plotdofs(3)),'}}$'],'interpreter','latex','FontSize',16);
else
    xlabel(varargin{1}{1},'interpreter','latex','FontSize',16);
    ylabel(varargin{1}{2},'interpreter','latex','FontSize',16);
    zlabel(varargin{1}{3},'interpreter','latex','FontSize',16);
end    

end
