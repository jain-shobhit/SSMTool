function figure_14_4
% Figure 14.4: Projecting the resonance surface shown in Fig. 14.3 onto the
% rho-epsilon plane results in an Arnol'd tongue (a). Combining the
% computation of resonance surfaces with the computation of quasi-periodic
% arcs, as shown in Section 9.2, enables the investigation of the so-called
% Arnol'd tongue scenario. This consists of a countable collection of
% Arnol'd tongues and a complementary Cantor-like set of quasi-periodic
% arcs. Panel (b) shows, from left to right, the 1:4, 3:11, 2:7, 3:10 and
% 1:3 tongues together with quasi-periodic arcs for rotation numbers
% obtained by a two-level recursive golden-mean subdivision of the interval
% [1/4,1/3]. The tongues are represented by a dot for each of the 1,500
% resonant orbits computed on each resonance surface. We observe that one
% of the quasi-periodic arcs enters the 2:7 tongue. This violates theory
% and is an artifact caused by the absence of control for discretization
% errors, a task undertaken in Part V.

% Note that the violation mentioned above is partially a result of
% round-off errors and thus machine-dependent.

if ~coco_exist('run1', 'run')
  run demo_resonant
end
if ~coco_exist('AT', 'run')
  run demo_arnold
end

figure(1)
clf
hold on
grid on
box on
view([0 90]);
axis([-0.12 0.67 -0.5 0.5]);

[tri X] = get_trisurf('run1');

C       = 0*ones(size(X,1), 1);
cmap    = repmat(linspace(0.55,0.95,100)', 1, 3);
colormap(cmap);
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'FaceColor', 'interp', ...
  'EdgeColor', 0.5*[1 1 1]);
surf([0 0;0 0], [0 0;0 0], [0 0;0 0], 'FaceColor', 'white', ...
  'EdgeColor', 'white', 'FaceAlpha', 0);

hold off

bd     = coco_bd_read({'AT' 'run_po'});
lab    = coco_bd_col(bd, 'LAB');
idx    = coco_bd_idxs(bd, 'AT');
lab_at = fliplr([lab{idx}]);
idx    = coco_bd_idxs(bd, 'QP');
lab_qp = fliplr([lab{idx}]);

figure(2)
clf
hold on
grid on
box on
view([0 90]);
axis([0.23 0.364 -0.1 0.1]);


surf([0 0;0 0], [0 0;0 0], [0 0;0 0], 'FaceColor', 'white', ...
  'EdgeColor', 'white', 'FaceAlpha', 0);

for lab=lab_qp
  runid = {'AT', sprintf('run_qp_%d',lab)};
  bd = coco_bd_read(runid);
  ro = coco_bd_col(bd, 'ro');
  ep = coco_bd_col(bd, 'eps');
  plot(ro, ep, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black');
end

for lab=lab_at
  runid = {'AT', sprintf('run_at_%d',lab)};
  bd = coco_bd_read(runid);
  ro = coco_bd_col(bd, 'ro');
  ep = coco_bd_col(bd, 'eps');
  plot(ro, ep, 'LineStyle', 'none', 'LineWidth', 2, ...
    'Color', [0.7 0.7 0.7], 'Marker', '.', 'MarkerSize', 9);
end

end

function [tri X] = get_trisurf(run)

atlas  = coco_bd_read(run, 'atlas');
charts = atlas.charts;
tri = [];
X   = [];
N   = numel(charts);
for k=1:N
  chart = charts{k};
  X     = [ X ; chart.x([153 154 1 2 3])' ]; %#ok<AGROW>
  ic    = [chart.nb chart.nb(1)];
  ix    = chart.id;
  for l=1:numel(ic)-1
    face = sort([ix ic(l) ic(l+1)]);
    if all(face>0) && (isempty(tri) ||  ~ismember(face, tri, 'rows'))
      tri = [tri ; face]; %#ok<AGROW>
    end
  end
end

end
