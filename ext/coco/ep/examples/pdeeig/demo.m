%% The Laplace operator
%
% We consider the boundary-value problem u_xx + u_yy = -la*(u + mu*exp(u))
% with u vanishing on the boundary of an L-shaped domain. We rely on a
% finite-difference discretization of the two-dimensional Laplacian and use
% continuation to compute eigenfunctions of the Laplacian on the given
% domain for mu = 0.

% Each figure shows the eigenfunction with normalized Euclidean norm
% ||u||=1.

%% Initial encoding

% The continuation problem encoded below includes two monitor functions
% that evaluate to the problem parameters mu and la, respectively, and four
% corresponding inactive continuation parameters 'mu' and 'la'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'la' is released and allowed to vary during
% continuation.

% Finite-difference discretization of the boundary-value problem defined
% using an inline, anonymous function. The order of discretization is
% parameterized by N. Specifically, the problem dimension equals
% (3N+1)*(2N+1)+const. It is worth experimenting with different values of N
% to explore convergence and computational efficiency.

N = 20;
P = 3*N+1;
Q = 2*N+1;
X = reshape(1:P*Q, P, Q);

Mask            = true(P,Q);
Mask(1:end,1)   = false;
Mask(1:end,end) = false;
Mask(1,1:end)   = false;
Mask(end,1:end) = false;
Mask(1:2*N+1,1:N+1) = false;

rows = X(Mask);
cols = rows;
o    = ones(numel(rows),1);
C    = sparse(rows, cols, o, P*Q, P*Q);
D    = sparse(rows, cols, -4*o, P*Q, P*Q);
cols = X(circshift(Mask, [0 1]));
D    = D + sparse(rows, cols, o, P*Q, P*Q);
cols = X(circshift(Mask, [0 -1]));
D    = D + sparse(rows, cols, o, P*Q, P*Q);
cols = X(circshift(Mask, [1 0]));
D    = D + sparse(rows, cols, o, P*Q, P*Q);
cols = X(circshift(Mask, [-1 0]));
D    = D + sparse(rows, cols, o, P*Q, P*Q);
D    = N^2*D;

rows = X(~Mask);
cols = rows;
o    = ones(numel(rows),1);
B    = sparse(rows, cols, o, P*Q, P*Q);

Id   = speye(P*Q, P*Q);

pdeeig    = @(u,p) D*u + p(2)*(C*u + p(1)*(C*exp(u))) - B*u;
pdeeig_dx = @(u,p) D + p(2)*(C + p(1)*(C*spdiags(exp(u),0,Id))) - B;
pdeeig_dp = @(u,p) [p(2)*(C*exp(u)) C*u+p(1)*(C*exp(u))];

p0 = [0; 0.1];
u0 = zeros(P*Q,1);

% Initialize continuation problem and settings associated with toolbox
% constructor and atlas algorithm.
prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', false);
% switch off bifurcation detection [memory limitation]
prob = coco_set(prob, 'ep', 'bifus', false);

%% Start continuation along trivial branch with mu=0 and u=0

% Branch points correspond to zero-amplitude eigenfunctions of the Laplace
% operator with corresponding eigenvalue -la.

prob = coco_set(prob, 'cont', 'PtMX', 200);

fprintf('\n Run=''%s'': Continue equilibria along trivial branch.\n', ...
  'pdeeig1');

bd = coco(prob, 'pdeeig1', 'ode', 'isol', 'ep', ...
  pdeeig, pdeeig_dx, pdeeig_dp, u0, {'mu' 'la'}, p0, ...
  1, 'la', [-5 50]);

%% Restart continuation along secondary branches of eigenfunctions

% The continuation problem encoded below includes two monitor functions
% that evaluate to the problem parameters mu and la, respectively, and four
% corresponding inactive continuation parameters 'mu' and 'la'. Its
% dimensional deficit equals 0. The call to the coco entry-point function
% indicates a desired manifold dimension of 1. To this end, the
% continuation parameter 'la' is released and allowed to vary during
% continuation.

labs = coco_bd_labs(bd, 'BP');

prob = coco_set(prob, 'cont', 'PtMX', [10 0]);
prob = coco_set(prob, 'cont', 'h0', 1);

for lab = labs
  % add toolbox outside call to coco
  prob2 = ode_BP2ep(prob, '', 'pdeeig1', lab);
  
  % add norm as regular monitor function
  [fdata, uidx] = coco_get_func_data(prob2, 'ep', 'data', 'uidx');
  xidx          = uidx(fdata.ep_eqn.x_idx);
  prob2         = coco_add_func(prob2, 'norm_x', @norm_x, [], ...
    'regular', 'norm_x', 'uidx', xidx);
  
  % add ||x||_2 = 1 as boundary event to compute unit eigenfunction
  prob2 = coco_add_event(prob2, 'UZ', 'BP', 'norm_x', 1);
  
  % compute eigenfunction
  runid = sprintf('pdeeig2_%02d', lab);
  fprintf(...
    '\n Run=''%s'': Switch branch at point %d in run ''%s''.\n', ...
    runid, lab, 'pdeeig1');

  bd2 = coco(prob2, runid, [], 1, {'la' 'norm_x'});
  
  % plot eigenfunction
  figure(1); clf
  coco_plot_sol(struct('plot_sol', @ef_plot, 'N', N), runid, '')
  view([-125 60])
  drawnow
end
