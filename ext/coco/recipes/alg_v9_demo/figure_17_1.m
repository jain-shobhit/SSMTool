function figure_17_1
% Figure 17.1: Covering of a manifold of equilibrium points of Example
% 17.1. The black dots mark events associated with psi_HB and the value 0,
% which are either Hopf points or neutral saddles. The circles mark fold
% points. From panel (a) it is evident that the loci of Hopf and fold
% points intersect each other transversally. When projected onto the
% (p_1,p_2) parameter plane, these curves have a tangency (b). The gray
% curve in panel (b) is the locus of zeros of psi_HB given by Eq. (17.4).
% The distinction between different events associated with the same monitor
% function is afforded by event handlers; see Fig. 17.3.

% Generate data
if ~(coco_exist('run', 'run'))
  run demo_alg_v9
end

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([0 0.5 0 0.25 0 10])
view([-15 15])

[atlas bd] = coco_bd_read('run', 'atlas', 'bd');
[tri X] = get_trisurf(atlas.charts, 3,4,1);

x0      = repmat([0 0 0], size(X,1), 1);
n       = [-0.5 -0.3 0.01];
C       = n*(X-x0)';
cmap    = repmat(linspace(0.55,0.95,100)', 1, 3);
colormap(cmap);
surf([0 0;0 0], [0 0;0 0], [0 0;0 0], 'FaceColor', 'white', ...
  'EdgeColor', 'white', 'FaceAlpha', 0)
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'FaceColor', 'interp', ...
  'EdgeColor', 0.5*[1 1 1], 'LineWidth', 0.5)

idx1 = coco_bd_idxs(bd, 'HB');
idx2 = coco_bd_idxs(bd, 'FO');
x  = coco_bd_col(bd, 'x');
p1 = coco_bd_col(bd, 'p1');
p2 = coco_bd_col(bd, 'p2');
plot3(p1(idx1)-0.001, p2(idx1)-0.001, x(idx1), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', '.', 'MarkerSize', 15, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black')
plot3(p1(idx2)-0.002, p2(idx2), x(idx2), 'LineStyle', 'none', ...
  'LineWidth', 2, 'Color', 'black', 'Marker', 'o', 'MarkerSize', 9.5, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
box on
axis([0 0.5 0 0.255])

x = linspace(0,0.25,100);
y = (1+sqrt(1-4*x)+2*x)./(4+2*x);
plot(y, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])
y = (1-sqrt(1-4*x)+2*x)./(4+2*x);
plot(y, x, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6])

plot(p1(idx1), p2(idx1), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 15, 'MarkerEdgeColor', ...
  'black', 'MarkerFaceColor', 'black')
plot(p1(idx2), p2(idx2), 'LineStyle', 'none', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', 'o', 'MarkerSize', 9.5, ...
  'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white')

hold off

end

function [tri X] = get_trisurf(charts, IX, IY, IZ)

tri = [];
X   = [];
N   = numel(charts);
for k=1:N
  chart = charts{k};
  X     = [ X ; chart.x([IX IY IZ])' ]; %#ok<AGROW>
  ic    = [chart.nb chart.nb(1)];
  ix    = chart.id;
  for l=1:numel(ic)-1
    face = sort([ix ic(l) ic(l+1)]);
    if all(face>0) && (isempty(tri) ||  ~ismember(face, tri, 'rows'))
      tri  = [tri ; face]; %#ok<AGROW>
    end
  end
end

end
