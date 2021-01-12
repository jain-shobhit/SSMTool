function figure_9_1
% Figure 9.1: Two patches of a lift onto R^2 of the fundamental domain of
% the transport equation in Eq. (9.13) next to each other with a sketch of
% the characteristic field (a). Applying the method of characteristics with
% initial condition u(theta_1,0) maps the closed curve u(.,0) onto the
% closed curve u(.,2*pi), whereby the parameterization experiences a
% rotation by 2*pi*rho. In state space, this translates to a flow on an
% invariant torus for which there exists a curve
% y(phi)=upsilon(phi,0):=u(2*pi*phi,0), for phi in [0,1], that is mapped
% onto itself under the original flow after time T=2*pi/omega_2 and under a
% rotation of the parameterization by rho, as illustrated in panel (b)

% Generate data
if ~coco_exist('run0', 'run')
  run demo_langford
end

coco_use_recipes_toolbox coll_v1 msbvp_v1

% Plot data: panel (b)
figure(1)
clf
hold on
grid on
axis([-1.5 1.5 -2 2 -0.5 2])
view([50 15])

plot_torus('run0', 1);

hold off

coco_use_recipes_toolbox

end

function [sol data] = plot_torus(run, lab, fac)

if nargin<3
  fac=1;
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
plot3(sol{30}.x(:,1), sol{30}.x(:,2), sol{30}.x(:,3), 'LineStyle', '-', ...
  'LineWidth', 2, 'Color', 'black');

end
