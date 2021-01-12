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
% quasiperiodic invariant tori. Selected members of this family are shown
% in panels (b) to (f).

%% Encoding

% The continuation problem structure encoded below includes three monitor
% functions that evaluate to the problem parameters, and two corresponding
% inactive continuation parameters 'om', 'c', and 'a'. Its dimensional
% deficit equals -1. The call to the coco entry-point function indicates a
% desired manifold dimension of 1. To this end, the two continuation
% parameters 'om' and 'a' are released and allowed to vary during
% continuation.

% Construct 'coll' arguments

% pnames = [  om     c     a ]
p0       = [ 1.5111; 0.11; 0.1 ];

T_po = 2*pi; % Approximate period
N    = 10;  % 2N+1 = Number of orbit segments
tout = linspace(0, T_po, 2*N+2);

T_ret = 2*pi/p0(1); % return time = 2*pi/om
tt    = linspace(0,1,10*(2*N+1))';
t1    = T_ret*tt;
coll  = cell(1,2*N+1);
for i=1:2*N+1 % For each point on orbit, flow for return time and reconstitute 3D trajectory
  x1       = [2*cos(tout(i)) 2*sin(tout(i))];
  [~, x1] = ode45(@(t,x) vdp(t,x,p0), 10*t1, x1); % Transient simulation
  [~, x1] = ode45(@(t,x) vdp(t,x,p0), t1, x1(end,:));
%    x1       = [2*cos(t1-tout(i)) -2*sin(t1-tout(i))];
  coll{i}  = {@vdp @vdp_DFDX @vdp_DFDP @vdp_DFDT t1 x1 p0};
end

% Construct boundary conditions data, Fourier transform and rotation matrix
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
%Fi = [ones(2*N+1,1)   reshape([cos(Th);sin(Th)], [2*N+1 2*N])];
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);

varrho = T_ret/T_po;
Th  = -(1:N)*2*pi*varrho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);

data = struct();
data.F  = kron(F, eye(2));
data.RF = kron(R*F, eye(2));

% Use the 'F+dF' option of the ode_isol2bvp constructor, since @torus_bc
% evaluates both the residual of the boundary conditions and the
% corresponding Jacobian. Turns off adaptivity to allow for graphing with
% surf function.
prob = coco_prob();
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'coll', 'NTST', 40);
prob = coco_set(prob, 'cont', 'NAdapt', 0, 'h_max', 2, 'PtMX', 200);
prob = coco_set(prob, 'corr', 'ItMX', 20);
prob = ode_isol2bvp(prob, '', coll, {'om' 'c' 'a'}, ...
  @torus_bc, data, 'F+dF');

fprintf('\n Run=''%s'': Continue family of quasiperiodic invariant tori.\n', ...
  'torus');

coco(prob, 'vdP_torus', [], 1, {'a' 'om'}, {[0.1 2] [1 2]});

%% Graphical representation of stored solutions

% Plot data: panel (a)
figure(1); clf
coco_plot_bd('vdP_torus', 'a', 'om')
grid on

% Plot data: panels (b)-(f)
labs = [1 4 6 8 11];

for i=1:5
  figure(i+1); clf; hold on; grid on
  
  [sol, data] = bvp_read_solution('', 'vdP_torus', labs(i));
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
