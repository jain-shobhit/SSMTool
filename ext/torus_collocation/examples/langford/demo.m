
%% Example 9.1 in recipes

% The continuation problem structure encoded below corresponds to a family
% of 151*(2N+1)+2+8N zero functions (120*(2N+1) collocation conditions,
% 27*(2N+1) continuity conditions, 2N+1 boundary conditions on the interval
% lengths, 3*(2N+1) Fourier interpolation conditions, two phase conditions,
% and 8N gluing conditions) in terms of 155*(2N+1) continuation variables
% (150*(2N+1) basepoint values, 2N+1 interval lengths, and 4*(2N+1) problem
% parameters), a family of four monitor functions, and the corresponding
% inactive continuation parameters 'om', 'ro', 'eps', 'T_ret'. Its
% dimensional deficit is -2. It follows that two continuation parameters
% must be released in order to obtain a unique solution to a closed
% continuation problem. Given 'om', and for 'eps'=0, a unique family of
% trajectory segments is obtained by solving also for 'ro', and 'T_ret'.
% Additional variations in 'eps' result in a one-dimensional branch of
% approximate quasiperiodic invariant tori. Alternatively, variations in
% 'ro' describe a one-dimensional branch of relative periodic orbits.

p0     = [3.5; 0.35; 0];
T_po   = 5.3; % Approximate period
N      = 50;  % 2N+1 = Number of orbit segments
tout   = linspace(0, T_po, 2*N+2);
[t x0] = ode45(@(t,x) lang_red(x,p0), tout, [0.3; 0.4]); % Approximate periodic orbit of reduced system

T_ret = 2*pi/p0(1); % return time = 2*pi/om
tt    = linspace(0,1,20*(2*N+1))';
t1    = T_ret*tt;
stt   = sin(tt*2*pi);
ctt   = cos(tt*2*pi);
xinit = zeros(numel(tt),3,2*N+1);
for i=1:2*N+1 % For each point on orbit, flow for return time and reconstitute 3D trajectory
  [t xx]    = ode45(@(t,x) lang_red(x,p0), [0 T_ret], x0(i,:));
  xx        = interp1(t, xx, t1);
  x1        = [ctt.*xx(:,1) stt.*xx(:,1) xx(:,2)]; % x0 -> [cos(t)*x0(:,1) sin(t)*x0(:,1) x0(:,2)]
  xinit(:,:,i) = x1; 
%   coll_args = [coll_args {@lang @lang_DFDX @lang_DFDP t1 x1 [p0; T_ret]}];
end
varrho = T_ret/T_po;

% construct continuation problem
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL', false);
prob = coco_set(prob, 'cont', 'NAdapt', 0, 'h_max', 10, 'PtMX', 40);
torargs = {@lang @lang_DFDX @lang_DFDP t1 xinit {'om','rho','eps','om1','om2','varrho'} [p0' 2*pi/T_po p0(1) varrho]};
% we choose negative om1 above because the rotation direction is opposite.
% In other words, the value of om1 can be both positive and negative
prob = ode_isol2tor(prob, '', torargs{:});

fprintf('\n Run=''%s'': Continue family of quasiperiodic invariant tori.\n', ...
  'torus');
coco(prob, 'run_eps', [], 1, {'eps','rho','om1','om2','om','varrho'}, [-0.3 0.3]);

%% new run from previous solution - varying varrho
bd  = coco_bd_read('run_eps');
lab = coco_bd_labs(bd, 'EP');
lab = max(lab);
prob = coco_prob();
prob = coco_set(prob, 'coll', 'MXCL', false);
prob = coco_set(prob, 'cont', 'NAdapt', 0, 'h_max', 10, 'PtMX', 50);
prob = ode_tor2tor(prob, '', 'run_eps', lab);
fprintf('\n Run=''%s'': Continue family of quasiperiodic invariant tori.\n', ...
  'run_oms');
coco(prob, 'run_oms', [], 1, {'om1','varrho','om','om2','eps','rho'}, [1 2]);

%% switch to secondary branch via branch-switching
bd  = coco_bd_read('run_oms');
lab = coco_bd_labs(bd, 'BP');
lab = min(lab);
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 0, 'h_max', 10, 'PtMX', 30);
prob = ode_BP2tor(prob, '', 'run_oms', lab);

fprintf('\n Run=''%s'': Continue family of quasiperiodic invariant tori.\n', ...
  'run_oms_BP');

coco(prob, 'run_oms_BP', [], 1, {'om1','varrho','om','om2','eps','rho'}, [1 2]);

%% plot results
% continuation path
figure;
coco_plot_bd('run_eps', 'rho', 'eps');hold on

figure; hold 
coco_plot_bd('run_oms', 'om1', 'varrho');hold on
coco_plot_bd('run_oms_BP', 'om1', 'varrho'); % secondary branch
grid on

% torus
plot_torus('','run_oms_BP', 1, [1 2 3]);
