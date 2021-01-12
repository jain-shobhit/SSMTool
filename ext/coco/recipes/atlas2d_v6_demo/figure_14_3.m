function figure_14_3
% Figure 14.3: A so-called resonance surface obtained by applying the
% 2-dimensional covering algorithm from Sect. 14.2 to a boundary-value
% problem for resonant periodic orbits of the Langford dynamical system,
% given by the vector field in Eq. (14.4). A projection of this surface
% onto the rho-epsilon-y_1(0) space is shown in panel (a). The striking
% cyclic triple-S shape is generic for 1:3 resonance surfaces. Similar
% shapes are observed for all m:n resonance surfaces. The apparent
% self-intersection is due to projection; the surface is locally isomorphic
% to a cylinder (b). A more familiar projection, a so-called Arnol'd
% tongue, is shown in Fig. 14.4(a).

% Generate data
if ~coco_exist('run1', 'run')
  run demo_resonant
end

[tri X] = get_trisurf('run1');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([-0.12 0.66 -0.5 0.5 0.24 1.45])
view([30 20])

x0      = repmat([0.3 0 0.8], size(X,1), 1);
n       = [1 -2 4];
C       = n*(X(:,[1 2 3])-x0)';
cmap    = repmat(linspace(0.55,0.95,100)', 1, 3);
colormap(cmap);
trisurf(tri, X(:,1), X(:,2), X(:,3), C, 'FaceColor', 'interp', ...
  'EdgeColor', 0.5*[1 1 1])
surf([0 0;0 0], [0 0;0 0], [0 0;0 0], 'FaceColor', 'white', ...
  'EdgeColor', 'white', 'FaceAlpha', 0)

hold off

% Plot data: panel (b)
figure(2)
clf
hold on
grid on
axis([-0.2 1.8 -0.5 0.5 0.24 1.45])
view([30 20])

x0      = repmat([0.3 0 0.8], size(X,1), 1);
n       = [2 -1 2];
C       = n*(X(:,[5 2 3])-x0)';
cmap    = repmat(linspace(0.55,0.95,100)', 1, 3);
colormap(cmap);
trisurf(tri, X(:,5), X(:,2), X(:,3), C, 'FaceColor', 'interp', ...
  'FaceAlpha', 0.8, 'EdgeColor', 0.5*[1 1 1])

hold off

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
