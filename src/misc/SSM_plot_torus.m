function SSM_plot_torus(qTr,dofs,varargin)
% SSM_PLOT_TORUS This function present visualization of a torus represented
% by a set of trajectories on such a torus. 
%
% SSM_PLOT_TORUS(QTR,DOFS,VARARGIN)
% VARARGIN = AXISLABELS
%
% qtr:        a set of trajectories stored in a 3D array with format [nt, dim, nsegs]
% dofs:       dofs for plotting in x-y-z directions
% axisLables: {xlabel,ylabel,zlabel} 
%
% See also: SSM_PO_READ_SOLUTION

assert(numel(dofs)==3, 'The size of dofs is not three. Plotting torus is disabled');
dim = size(qTr,2);
assert(dim>=max(dofs), 'The dimension of trajectory is not enohgh or dofs are too large');
ytube = permute(qTr,[3 1 2]);
figure; hold on
h = surf(ytube(:,:,dofs(1)),ytube(:,:,dofs(2)),ytube(:,:,dofs(3)));
set(h,'edgecolor','none')
set(h,'FaceAlpha',0.5);
view([1,1,1])
grid on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
if isempty(varargin)
    xlabel(['$x_{\mathrm{',num2str(dofs(1)),'}}$'],'interpreter','latex','FontSize',16);
    ylabel(['$x_{\mathrm{',num2str(dofs(2)),'}}$'],'interpreter','latex','FontSize',16);
    zlabel(['$x_{\mathrm{',num2str(dofs(3)),'}}$'],'interpreter','latex','FontSize',16);
else
    xlabel(varargin{1}{1},'interpreter','latex','FontSize',16);
    ylabel(varargin{1}{2},'interpreter','latex','FontSize',16);
    zlabel(varargin{1}{3},'interpreter','latex','FontSize',16);
end    

end