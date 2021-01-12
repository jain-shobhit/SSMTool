%% Continuation of invariant circle
%
% Continuation problem is equivalent to analysis of quasiperiodic invariant
% torus for vector field in torus.m. Here, an explicit solution is used to
% eliminate the need to invoke multiple instances of the 'coll' toolbox.

%% Encoding

% The continuation problem structure encoded below includes two monitor
% functions that evaluate to the problem parameters, and two corresponding
% inactive continuation parameters 'Om' and 'om'. Its dimensional deficit
% equals -1. The call to the coco entry-point function indicates a desired
% manifold dimension of 1. To this end, both continuation parameters are
% released and allowed to vary during continuation.

% Construct family of algebraic zero problems mapping initial points to
% final points
om   = 1.5;
Om   = 1;
N    = 15;  % 2N+1 = Number of points on invariant circle
vphi = 2*pi*linspace(0,1,2*N+2);
T    = 2*pi/om;
prob = coco_prob();
z0_idx = zeros(2*(2*N+1),1);
z1_idx = zeros(2*(2*N+1),1);
T_idx  = zeros(2*N+1,1);
for i=1:2*N+1
  off = 7*(i-1);
  z0 = repmat((1+om^2)/om^2, [2 1]).*[cos(vphi(i)); sin(vphi(i))];
  z1 = repmat((1+om^2)/om^2, [2 1]).*[cos(T+vphi(i)); sin(T+vphi(i))];
  prob = coco_add_func(prob, sprintf('f%d', i), @map, [], 'zero', ...
    'u0', [z0; z1; om; Om; 2*pi/om]);
  z0_idx(2*(i-1)+(1:2)') = 7*(i-1)+(1:2)';
  z1_idx(2*(i-1)+(1:2)') = 7*(i-1)+(3:4)';
  T_idx(i) = 7*i;
end

% Glue parameters
for i=1:2*N
  prob = coco_add_glue(prob, sprintf('glue%d',i), 7*i+(5:6), 5:6);
end
prob = coco_add_pars(prob, 'pars', [5;6], {'om' 'Om'});

% Construct boundary conditions data, Fourier transform and rotation matrix
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1); % 2N+1 by 2N+1

varrho = 1/1.5111;
Th  = (1:N)*2*pi*varrho;
SIN = [ zeros(size(Th)) ; sin(Th) ];
R   = diag([1 kron(cos(Th), [1 1])]);
R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);

data = struct();
data.F  = kron(F, eye(2));
data.RF = kron(R*F, eye(2));

% Append boundary conditions
data.T_idx  = 1:(2*N+1);
data.x0_idx = 2*N+1+(1:2*(2*N+1));
data.x1_idx = 6*N+3+(1:2*(2*N+1));
data.p_idx  = 10*N+5+(1:2);
prob = coco_add_func(prob, 'bc', @bc, data, 'zero', ...
  'uidx', [T_idx; z0_idx; z1_idx; 5;6], 'F+dF');
prob = coco_set(prob, 'cont', 'h_max', 10);

fprintf('\n Run=''%s'': Continue family of invariant circles.\n', ...
  'inv_circ');

coco(prob, 'inv_circ', [], 1, {'om' 'Om'}, [0.5 1.5]);

%% Graphical representation of stored solutions

figure(1); clf
coco_plot_bd('inv_circ', 'om', 'Om')
grid on
