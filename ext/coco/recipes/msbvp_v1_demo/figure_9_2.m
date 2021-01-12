function figure_9_2
% Figure 9.2: A continuation of quasi-periodic invariant tori of the
% Langford system in Eq. (9.4) results in the curve shown in panel (a),
% which is referred to as a quasi-periodic arc or a quasi-periodic hair.
% Selected members of this family are shown in panels (b) to (f). The
% accumulation of orbits on the torus suggests that this family approaches
% the vicinity of a 1:3 resonance at both ends of the arc. The labels
% correspond to the session output included in the text.

% Generate data
if ~coco_exist('run_eps', 'run')
  run demo_langford
end

coco_use_recipes_toolbox coll_v1 msbvp_v1

% Extract data
bd  = coco_bd_read('run_eps');
eps = coco_bd_col(bd, 'eps');
ro  = coco_bd_col(bd, 'ro');

% Plot data: panel (a)
figure(1)
clf
hold on
grid on
axis([0.37 0.4 -0.05 0.05])

plot(ro, eps, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'black', ...
  'Marker', '.', 'MarkerSize', 12)

hold off

% Plot data: panels (b)-(f)
labs = [6, 4, 1, 10, 12];

for i=1:5
figure(i+1)
clf
hold on
grid on
axis([-1.5 1.5 -2 2 -0.5 2])
view([50 15])

plot_torus('run_eps', labs(i))

hold off
end 

coco_use_recipes_toolbox

end

function [sol data] = plot_torus(run, lab, fac)

if nargin<3
  fac=0.75;
end

[sol data] = msbvp_read_solution('', run, lab);
N  = data.nsegs;
M  = ceil(fac*size(sol{1}.x,1));
x0 = zeros(N+1,3);
x1 = zeros(N+1,3);
XX = zeros(M,N+1);
YY = XX;
ZZ = XX;
for i=1:N+1
  n       = mod(i-1,N)+1;
  XX(:,i) = sol{n}.x(1:M,1);
  YY(:,i) = sol{n}.x(1:M,2);
  ZZ(:,i) = sol{n}.x(1:M,3);
  x0(i,:) = sol{n}.x(1,:);
  x1(i,:) = sol{n}.x(M,:);
end
surf(XX, YY, ZZ, 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.7, ...
  'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
  'LineWidth', 0.5);
plot3(x0(:,1), x0(:,2), x0(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
plot3(x1(:,1), x1(:,2), x1(:,3), 'LineStyle', '-', 'LineWidth', 2, ...
  'Color', 'black', 'Marker', '.', 'MarkerSize', 12);

end
