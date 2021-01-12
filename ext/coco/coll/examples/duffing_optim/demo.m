%% Stationary points in Duffing oscillator
%
% We apply the method of staged continuation to the search for stationary
% values of amp along a family of periodic solutions of the harmonically
% excited Duffing oscillator
%
%     x'' + 2*zeta*x' + x + alpha*x^3 = amp*cos(omega*t)
%
% In the first stage of continuation, a local extremum in amp is detected
% as a fold along the solution manifold. By construction, this coincides
% with a branch point of the augmented continuation problem that includes
% the adjoint conditions associated with the zero and monitor functions. In
% the second stage, we continue along the secondary branch until the
% Lagrange multiplier associated with amp equals 1. As described in the
% tutorial document, at this point, 'amp' = 2*'alpha'*'d.alpha'

%% initial solution guess
p0 = [0.1; 0.5; 0.06; 1.2];
t0 = linspace(0,2*pi/p0(4),500)';
x0 = [p0(3)/(1-p0(4)^2)*cos(p0(4)*t0) ...
  -p0(3)*p0(4)/(1-p0(4)^2)*sin(p0(4)*t0) ...
  p0(4)*t0];

prob = coco_prob;
prob = coco_set(prob, 'cont', 'h_max', 10);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'NPR', 10);

%% first run to find initial fold

prob1 = coco_set(prob, 'coll', 'NTST', 30);

% zero problems
coll_args = {@duffing, @duffing_dx, @duffing_dp, @duffing_dxdx,...
  @duffing_dxdp, @duffing_dpdp, t0, x0, p0};
bvp_args = [coll_args, {{'zeta' 'alpha' 'amp' 'omega'}, @duffing_bc, ...
  @duffing_bc_du, @duffing_bc_dudu}];
prob1 = ode_isol2bvp(prob1, '', bvp_args{:});

% adjoints
prob1 = adjt_isol2bvp(prob1, '');

% continuation
cont_pars = {'amp' 'd.amp' 'd.zeta' 'd.alpha' 'd.omega'};
coco(prob1, 'duffing1', [], 1, cont_pars, [0.06 0.5]);
 
%% switch at fold to branch with nontrivial multipliers
bd1 = coco_bd_read('duffing1');
BPlab = coco_bd_labs(bd1, 'BP');

% zero problems
prob2 = ode_BP2bvp(prob, '', 'duffing1', BPlab(1));    

% adjoints
prob2 = adjt_BP2bvp(prob2, '', 'duffing1', BPlab(1));

% continuation
cont_pars = {'d.amp' 'amp' 'd.zeta' 'd.alpha' 'd.omega'};
coco(prob2, 'duffing2', [], 1, cont_pars, {[0 1] [0.1 2]});
