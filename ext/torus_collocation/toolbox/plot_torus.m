function plot_torus(toid, run, lab, outdof, varargin)
% PLOT_TORUS This function plots torus characterized by outdof
%
% PLOT_TORUS(TOID, RUN, LAB, OUTDOF, VARARGIN)
% VARARGIN = [FAC], [AXLABS]
%
% toid   - the instance of torus, a typical used one is ''
% run    - the id of run
% lab    - the label of solution
% outdof - dofs for torus plotting. If the number of dofs is two, the plot
%          of dof1-time-dof2 will be given. If the number of dofs is three,
%          the plot of dof1-dof2-dof3 will be given
% fac    - fraction of torus plot
% axLabs - {xlabel,zlabel} or {xlabel,ylabel,zlabel}
%
% See also: TOR_READ_SOLUTION

labels = [];
if nargin<5
    fac= 1;
else
    if isnumeric(varargin{1})
        fac = varargin{1};
        if nargin==6
            labels = varargin{2};
        end
    else
        fac = 1;
        labels = varargin{1};
    end
end

ndof = numel(outdof);
assert(ndof<4, 'visualization of 4 or higher dimensional torus is not supported');
assert(ndof>1, 'the dimension of outdof should be larger than one');

[sol, data] = tor_read_solution(toid, run, lab);
Tube = sol.xbp; % [nt dim nsegs+1]
nt   = size(Tube,1);
M    = ceil(nt*fac);
Tube = permute(Tube, [1,3,2]);
Tube = Tube(1:M,:,:); % a fraction of torus
Tube = Tube(:,1:end-1,outdof);
Nseg = size(Tube,2);
% Fourier expansion to recover
nsegs = 129; % 128+1
angs  = linspace(0,2*pi,nsegs);
tube  = Fourier_recover(permute(Tube,[2,1,3]), data.bc_data.Fs, angs);
tube  = permute(tube, [2,1,3]);
if ndof<3
    % plot of x1-t-x2
    % recovered tube
    tbp = repmat(sol.tbp, [1,nsegs]);
    tbp = tbp(1:M,:);
    figure; hold on
    surf(tube(:,:,1), tbp, tube(:,:,2), 'FaceAlpha', 0.7, ... %, 'FaceColor', 0.9*[1 1 1]
        'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 'none', ...
        'LineWidth', 0.5); % 
    plot3(tube(1,:,1), zeros(nsegs,1), tube(1,:,2), 'LineStyle', '-', 'LineWidth', 2, ...
        'Color', 'black');
    plot3(tube(end,:,1), tbp(end)*ones(nsegs,1), tube(end,:,2), 'LineStyle', '-', 'LineWidth', 2, ...
        'Color', 'blue');
    % trajectories at computed angles
    tbp = repmat(sol.tbp, [1,Nseg]);
    tbp = tbp(1:M,:);    
    for i=1:Nseg
        plot3(Tube(:,i,1), tbp, Tube(:,i,2), 'LineStyle', '-', ...
            'Color', 0.6*[1 1 1], 'LineWidth', 0.5);
    end
    plot3(Tube(1,:,1), zeros(Nseg,1), Tube(1,:,2), 'LineWidth', 2, ...
        'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
    plot3(Tube(end,:,1), tbp(end)*ones(Nseg,1), Tube(end,:,2), 'LineWidth', 2, ...
        'Color', 'blue', 'Marker', '.', 'MarkerSize', 12);    
    view([50 15]); grid on
    ylabel('$t$','interpreter','latex','FontSize',14);
    if isempty(labels)
        xlab = ['$x_\mathrm{',num2str(outdof(1)),'}$'];
        zlab = ['$x_\mathrm{',num2str(outdof(2)),'}$'];
    else
        xlab = labels{1};
        zlab = labels{2};
    end
    xlabel(xlab,'interpreter','latex','FontSize',14);
    zlabel(zlab,'interpreter','latex','FontSize',14);
    title(['$\omega_1=',num2str(sol.p(end-2)),', \omega_2=',num2str(sol.p(end-1)),'$'],...
        'interpreter','latex','FontSize',14);
    hold off
else
    % plot of x1-x2-x3
    % recovered tube
    figure; hold on
    surf(tube(:,:,1), tube(:,:,2), tube(:,:,3), 'FaceAlpha', 0.7, ... %, 'FaceColor', 0.9*[1 1 1]
        'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 'none', ...
        'LineWidth', 0.5);
    plot3(tube(1,:,1), tube(1,:,2), tube(1,:,3), 'LineStyle', '-', 'LineWidth', 2, ...
        'Color', 'black');
    plot3(tube(end,:,1), tube(end,:,2), tube(end,:,3), 'LineStyle', '-', 'LineWidth', 2, ...
        'Color', 'blue');
    % trajectories at computed angles
    for i=1:Nseg
        plot3(Tube(:,i,1), Tube(:,i,2), Tube(:,i,3), 'LineStyle', '-', ...
            'Color', 0.6*[1 1 1], 'LineWidth', 0.5);
    end
    plot3(Tube(1,:,1), Tube(1,:,2), Tube(1,:,3), 'LineWidth', 2, ...
        'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
    plot3(Tube(end,:,1), Tube(end,:,2), Tube(end,:,3), 'LineWidth', 2, ...
        'Color', 'blue', 'Marker', '.', 'MarkerSize', 12);        
    view([50 15]); grid on
    if isempty(labels)
        xlab = ['$x_\mathrm{',num2str(outdof(1)),'}$'];
        ylab = ['$x_\mathrm{',num2str(outdof(2)),'}$'];
        zlab = ['$x_\mathrm{',num2str(outdof(3)),'}$'];
    else
        xlab = labels{1};
        ylab = labels{2};
        zlab = labels{3};
    end
    xlabel(xlab,'interpreter','latex','FontSize',14);
    ylabel(ylab,'interpreter','latex','FontSize',14);
    zlabel(zlab,'interpreter','latex','FontSize',14);
    title(['$\omega_1=',num2str(sol.p(end-2)),', \omega_2=',num2str(sol.p(end-1)),'$'],...
        'interpreter','latex','FontSize',14);
    hold off
end

end
