%% A Brusselator model
%
% We consider a simplified Brusselator model consisting of the two
% diffusively coupled ordinary differential equations
%
%   u_t = delta*u_xx + alpha + u^2*v - (beta + 1)*u,
%   v_t = rho*delta*v_xx + beta*u - u^2*v
%
% and the associated boundary conditions
%
%    u(t,0) = u(t,1) = alpha, v(t,0) = v(t,1) = beta/alpha
%
% In looking for equilibrium solutions, we rely on a finite-difference
% discretization of the one-dimensional Laplacian and perform several
% continuation runs of equilibria and associated bifurcations.

% Figure 1 shows a one-dimensional family of equilibria under variations in
% beta, for fixed alpha, delta, and rho. Figure 2 shows one-dimensional
% families of saddle-node and hopf equilibria under simultaneous variation
% in beta and delta.

%% Initial encoding

% The continuation problem encoded below includes four monitor functions
% that evaluate to the problem parameters alpha, beta, delta, and rho,
% respectively, and four corresponding inactive continuation parameters
% 'al', 'be', 'de', and 'ro'. Its dimensional deficit equals 0. The call to
% the coco entry-point function indicates a desired manifold dimension of
% 1. To this end, the continuation parameter 'be' is released and allowed
% to vary during continuation.

% Finite-difference discretization of the Brusselator equilibrium
% boundary-value problem defined using an inline, anonymous function. The
% order of discretization is parameterized by N.

N  = 20;
X  = 1:N+1;
Y  = N+1+X;
B  = [1; zeros(N-1,1); 1];
BX = repmat(B,1,N+1);
BP = repmat(B,1,4);
C  = [0; ones(N-1,1); 0];
CX = repmat(C,1,N+1);
CP = repmat(C,1,4);
D  = diag([0 -2*ones(1,N-1) 0]) + ...
  diag([ones(1,N-1) 0],-1) + ...
  diag([0 ones(1,N-1)],1);
D  = N^2*D;
ID = eye(N+1,N+1);
O  = ones(N+1,1);
ZE = zeros(N+1,1);

bruss    = @(u,p) [
  C.*(p(3)*D*u(X) + p(1) + u(X).^2.*u(Y) - (p(2)+1)*u(X)) + B.*(p(1)-u(X))
  C.*((p(3)*p(4))*D*u(Y) + p(2)*u(X) - u(X).^2.*u(Y)) + B.*(p(2)/p(1)-u(Y))
  ];
bruss_dx = @(u,p) [
  CX.*(p(3)*D + 2*diag(u(X).*u(Y)) - (p(2)+1)*ID) - BX.*(ID) ...
  CX.*(diag(u(X).^2))
  CX.*(p(2)*ID - 2*diag(u(X).*u(Y))) ...
  CX.*((p(3)*p(4))*D - diag(u(X).^2)) - BX.*(ID)
  ];
bruss_dp = @(u,p) [
  C+B             C.*(-u(X))       C.*(D*u(X))                ZE
  -B.*(p(2)/p(1)^2)  C.*(u(X))+B.*(1/p(1))  C.*(p(4)*D*u(Y))  C.*(p(3)*D*u(Y))
  ];

p0 = [1; 3; 0.075; 1];
u0 = [p0(1)*ones(N+1,1); (p0(2)/p0(1))*ones(N+1,1)];

% Initialize continuation problem and settings associated with toolbox
% constructor and atlas algorithm.

prob = coco_prob();
prob = coco_set(prob, 'ode', 'vectorized', false);

%% Start continuation of equilibria from constant initial solution guess

% Include toolbox-specific 'ep.test.SN', 'ep.test.HB', 'ep.test.USTAB', and
% 'atlas.test.FP' continuation parameters for visual inspection of
% saddle-node, hopf, instability, and fold point monitor functions.

fprintf('\n Run=''%s'': Continue equilibria along primary branch.\n', ...
  'brus1');

bd1 = coco(prob, 'brus1', 'ode', 'isol', 'ep', ...
  bruss, bruss_dx, bruss_dp, u0, {'al' 'be' 'de' 'ro'}, p0,  ...               % ep toolbox arguments
  1, {'be' 'ep.test.SN' 'ep.test.HB' 'ep.test.USTAB' 'atlas.test.FP'}, [2 7]); % cont toolbox arguments

%% Restart continuation of equilbria along secondary branches

BPlabs = coco_bd_labs(bd1, 'BP');

for lab = BPlabs
  
  runid = sprintf('brus2_%02d', lab);
  fprintf('\n Run=''%s'': Switch branch at point %d in run ''%s''.\n', ...
    runid, lab, 'brus1');

  coco(prob, runid, 'ode', 'BP', 'ep', ...
    'brus1', lab, ...                                                            % ep toolbox arguments
    1, {'be' 'ep.test.SN' 'ep.test.HB' 'ep.test.USTAB' 'atlas.test.FP'}, [2 7]); % cont toolbox arguments
end

%% Start continuation along solutions to saddle-node and hopf equilibrium problem

% The continuation problems encoded below includes four monitor functions
% that evaluate to the problem parameters alpha, beta, delta, and rho,
% respectively, and four corresponding inactive continuation parameters
% 'al', 'be', 'de', and 'ro'. Their dimensional deficits equal -1. The
% calls to the coco entry-point function indicates a desired manifold
% dimension of 1. To this end, the continuation parameters 'de' and 'be'
% are both released and allowed to vary during continuation.

% Start from SN point on secondary branch

rrun   = sprintf('brus2_%02d', BPlabs(1));
bd     = coco_bd_read(rrun);
SNlabs = coco_bd_labs(bd, 'SN');
vals   = coco_bd_vals(bd, SNlabs, 'be');
[v,i]  = min(vals); %#ok<ASGLU>
rlab   = SNlabs(i);

fprintf(...
  '\n Run=''%s'': Continue SN points from point %d in run ''%s''.\n', ...
  'brus_SN', rlab, rrun);

prob = coco_set(prob, 'cont', 'PtMX', 200);
coco(prob, 'brus_SN', 'ode', 'SN', 'SN', ...
  rrun, rlab, ...                    %  ep toolbox arguments
  1, {'de' 'be' }, {[0 0.2] [2 7]}); % cont toolbox arguments

% Start from HB at primary branch
HBlabs = coco_bd_labs(bd1, 'HB');
rlab   = HBlabs(1);

fprintf(...
  '\n Run=''%s'': Continue HB points from point %d in run ''%s''.\n', ...
  'brus_HB1', rlab, 'brus1');

coco(prob, 'brus_HB1', 'ode', 'HB', 'HB', ...
  'brus1', rlab, ...                             %  ep toolbox arguments
  1, {'de' 'be' 'ep.test.BT'}, {[0 0.2] [2 7]}); % cont toolbox arguments

% Start from HB at secondary branch
rrun   = sprintf('brus2_%02d', BPlabs(1));
bd     = coco_bd_read(rrun);
HBlabs = coco_bd_labs(bd, 'HB');
vals   = coco_bd_vals(bd, HBlabs, 'be');
[v,i]  = min(vals);
rlab   = HBlabs(i);

fprintf(...
  '\n Run=''%s'': Continue HB points from point %d in run ''%s''.\n', ...
  'brus_HB2', rlab, rrun);

coco(prob, 'brus_HB2', 'ode', 'HB', 'HB', ...
  rrun, rlab, ...                    %  ep toolbox arguments
  1, {'de' 'be' }, {[0 0.2] [2 7]}); % cont toolbox arguments

%% Graphical representation of stored solutions

% "One-parameter continuation" of equilibria
figure(1); clf; hold on
thm = struct('special', {{'BP' 'HB'}});
coco_plot_bd(thm, 'brus1')
thm = struct('special', {{'FP' 'HB'}});
for lab=BPlabs
  coco_plot_bd(thm, sprintf('brus2_%02d', lab))
end
hold off; grid on

% "Two-parameter continuation" of saddle-node and Hopf bifurcations
figure(2); clf; hold on
coco_plot_bd('brus_SN', 'be', 'de', '||x||_2')
thm = struct('special', {{'BTP' 'BP'}});
coco_plot_bd(thm, 'brus_HB1', 'be', 'de', '||x||_2')
coco_plot_bd(thm, 'brus_HB2', 'be', 'de', '||x||_2')
hold off; grid on; view(3)
