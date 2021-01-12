%% Continuation of quasiperiodic invariant tori
%
% Two-dimensional quasiperiodic invariant tori are characterized by a
% parallel flow with irrational rotation number. In terms of a suitably
% defined torus function, the dynamics on the torus may be described by an
% invariant circle that is mapped to itself by the flow, such that the flow
% is equivalent to a rigid rotation on the circle.
%
% We obtain an approximate description of a quasiperiodic invariant torus
% in terms of a Fourier representation of the invariant circle and a finite
% collection of trajectory segments based at points on the invariant
% circle. The rigid rotation imposes an all-to-all coupled system of
% boundary conditions on the collection of trajectory segments.

% Panel (a) shows a quasiperiodic arc, representing a family of
% quasiperiodic invariant tori. The slope of the arc is the rotation number
% which is held constant during continuation. Selected members of this
% family are shown in panels (b) to (f).

%% Encoding

% The continuation problem structure encoded below includes two monitor
% functions that evaluate to the problem parameters, and two corresponding
% inactive continuation parameters 'Om' and 'om'. Its dimensional deficit
% equals -1. The call to the coco entry-point function indicates a desired
% manifold dimension of 1. To this end, both continuation parameters are
% released and allowed to vary during continuation.

% Construct 'coll' arguments
om   = 1.5;
Om   = 1;
N    = 15;  % 2N+1 = Number of orbit segments
vphi = 2*pi*linspace(0,1,2*N+2);
tau  = 2*pi/om*linspace(0,1,10*(2*N+1))';
rho  = (1+om^2)./(1+om^2-cos(om*tau)-om*sin(om*tau));
coll = cell(1,2*N+1);
for i=1:2*N+1
  up = repmat(rho, [1 2]).*[cos(Om*tau+vphi(i)) sin(Om*tau+vphi(i))];
  coll{i} = {@torus @torus_DFDX @torus_DFDP @torus_DFDT tau up [om Om]};
end

% Construct boundary conditions data, Fourier transform and rotation matrix
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);

varrho = 1/1.51111;
Th  = (1:N)*2*pi*varrho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);

data    = struct();
data.F  = kron(F, eye(2));
data.RF = kron(R*F, eye(2));

% Use the 'F+dF' option of the ode_isol2bvp constructor, since @torus_bc
% evaluates both the residual of the boundary conditions and the
% corresponding Jacobian.
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = coco_set(prob, 'cont', 'NAdapt', 1, 'h_max', 10);
prob = ode_isol2bvp(prob, '', coll, {'om' 'Om'}, @torus_bc, data, 'F+dF');

fprintf('\n Run=''%s'': Continue family of quasiperiodic invariant tori.\n', ...
  'torus');

coco(prob, 'torus', [], 1, {'om' 'Om'}, [0.5 1.5]);

%% Graphical representation of stored solutions

% Figure 1
figure(1); clf
coco_plot_bd('torus', 'om', 'Om')
grid on

% Plot data: panels (b)-(f)
bd = coco_bd_read('torus');
labs = coco_bd_labs(bd);

for i=1:numel(labs)
  figure(i+1); clf; hold on; grid on
  
  [sol, data] = bvp_read_solution('', 'torus', labs(i));
  N  = data.nsegs;
  M  = size(sol{1}.xbp,1);
  x0 = zeros(N+1,2);
  x1 = zeros(N+1,2);
  XX = zeros(M,N+1);
  YY = XX;
  ZZ = XX;
  for j=1:N+1
    n       = mod(j-1,N)+1;
    XX(:,j) = sol{n}.xbp(1:M,1);
    ZZ(:,j) = sol{n}.xbp(1:M,2);
    YY(:,j) = sol{n}.tbp(1:M);
    x0(j,:) = sol{n}.xbp(1,:);
    x1(j,:) = sol{n}.xbp(M,:);
  end
  surf(XX, YY, ZZ, 'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.7, ...
    'MeshStyle', 'column', 'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
    'LineWidth', 0.5);
  plot3(x0(:,1), zeros(N+1,1), x0(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  plot3(x1(:,1), sol{n}.T*ones(N+1,1), x1(:,2), 'LineStyle', '-', 'LineWidth', 2, ...
    'Color', 'black', 'Marker', '.', 'MarkerSize', 12);
  
  hold off; view([50 15])
end
