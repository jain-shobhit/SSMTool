% This demo proceeds in three steps:
%   1. Start at Hopf point and continue until homoclinic orbit is close.
%   2. Perform surgery to obtain initial guess for high-period orbit.
%   2. Iterate the remeshing map a number of times.
% The demo is executed three times, once without adaptation, and once for
% each moving-mesh algorithm.

%% Run with uniform mesh

coco_use_recipes_toolbox po_v1 coll_v3 % Add the po_v1 and coll_v3 toolboxes to the search path

t0 = (0:2*pi/100:2*pi)';
x0 = 0.01*(cos(t0)*[1 0 -1]-sin(t0)*[0 1 0]);
p0 = [0;6];

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 50);
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);
coco(prob, 'run1', [], 1, {'p1' 'po.period' 'po.seg.coll.NTST' 'po.seg.coll.err'}, [-1 1]);

sol1 = po_read_solution('', 'run1', 6); % Periodic orbit with T=1.8551e+01
f = marsden(sol1.x', repmat(sol1.p, [1 size(sol1.x,1)])); % Evaluate vector field at basepoints
[mn idx] = min(sqrt(sum(f.*f,1))); % Find basepoint closest to equilibrium

% Perform surgery on periodic orbit
scale = 25;
T     = sol1.t(end);
t0 = [sol1.t(1:idx,1) ; T*(scale-1)+sol1.t(idx+1:end,1)]; % Crank up period by factor scale
x0 = sol1.x;
p0 = sol1.p;

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', ceil(scale*50)); % Increase mesh resolution by factor scale
prob = coco_set(prob, 'coll', 'TOL', 1.0e-1);
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);

prob = coco_xchg_pars(prob, 'p2', 'po.period'); % Exchange parameters to fix period and free p2
coco(prob, 'run2', [], 0, {'p1' 'p2' 'po.period' 'po.seg.coll.NTST' 'po.seg.coll.err'});

%% Run with fixed-order moving mesh algorithm

coco_use_recipes_toolbox po_v3 coll_v5 % Add the po_v3 and coll_v5 toolboxes to the search path

t0 = (0:2*pi/100:2*pi)';
x0 = 0.01*(cos(t0)*[1 0 -1]-sin(t0)*[0 1 0]);
p0 = [0;6];

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', 20);
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);
prob = coco_set(prob, 'cont', 'NAdapt', 1); % Number of remesh-correct cycles
coco(prob, 'run1a', [], 1, {'p1' 'po.period' 'po.seg.coll.NTST' 'po.seg.coll.err'}, [-1 1]);

sol1 = po_read_solution('', 'run1a', 6); % Periodic orbit with T=1.9299e+01
f = marsden(sol1.x', repmat(sol1.p, [1 size(sol1.x,1)])); % Evaluate vector field at basepoints
[mn idx] = min(sqrt(sum(f.*f,1))); % Find basepoint closest to equilibrium

% Perform surgery on periodic orbit: note that reproduction of initial mesh
% in 'coll' constructor is exploited here to obtain a reasonable first
% solution point, our surgery inserts exactly NTSTINC mesh intervals such
% that the mesh of the excursion from the equilibrium remains unchanged
% (interpolation would introduce large errors), note also that the new
% period obtained here is significantly larger than above.

fac     = 50;
NTSTINC = 100;
NTST    = (numel(sol1.t)-1)/4; % NCOL = 4.
PTINC   = NTSTINC*4+2; % NCOL = 4.
T       = sol1.t(end);
scale   = fac*NTSTINC;
ti = linspace(sol1.t(idx,1), T*(scale-1)+sol1.t(idx+1,1), PTINC);
xi = zeros(PTINC,3);
for i=1:3
  xi(:,i) = linspace(sol1.x(idx,i), sol1.x(idx+1,i), PTINC);
end
t0 = [sol1.t(1:idx-1,1) ; ti' ; T*(scale-1)+sol1.t(idx+2:end,1)];
x0 = [sol1.x(1:idx-1,:) ; xi ; sol1.x(idx+2:end,:)];
p0 = sol1.p;

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', NTST+NTSTINC);
% Important: set TOL to something high --> MXCL does not occur, while we
% obtain good mesh
prob = coco_set(prob, 'coll', 'TOL', 1.0e-1);
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);

prob = coco_xchg_pars(prob, 'p2', 'po.period'); % Exchange parameters to fix period and free p2
prob = coco_set(prob, 'cont', 'NAdapt', 100); % Number of remesh-correct cycles
coco(prob, 'run2a', [], 0, {'p1' 'p2' 'po.period' 'po.seg.coll.NTST' 'po.seg.coll.err'}); % Iterate adaptation map

% for plot of dynamics of adaptation map use
%   bd = coco_bd_read('run2a');
%   plot(coco_bd_col(bd,'PT'), coco_bd_col(bd,'po.seg.coll.err'))

%% Run with variable-order moving mesh algorithm

coco_use_recipes_toolbox po_v3 coll_v6 % Add the po_v3 and coll_v6 toolboxes to the search path

t0 = (0:2*pi/100:2*pi)';
x0 = 0.01*(cos(t0)*[1 0 -1]-sin(t0)*[0 1 0]);
p0 = [0;6];

prob = coco_prob();
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);
prob = coco_set(prob, 'cont', 'NAdapt', 1); % Number of remesh-correct cycles
coco(prob, 'run1b', [], 1, {'p1' 'po.period' 'po.seg.coll.NTST' 'po.seg.coll.err'}, [-1 1]);

% find mesh point closest to equilibrium
sol1 = po_read_solution('', 'run1b', 6); % Periodic orbit with T=1.9277e+01
f = marsden(sol1.x', repmat(sol1.p, [1 size(sol1.x,1)])); % Evaluate vector field at basepoints
[mn idx] = min(sqrt(sum(f.*f,1))); % Find basepoint closest to equilibrium

% Perform surgery on periodic orbit: note that reproduction of initial mesh
% in coll constructor is exploited here to obtain a reasonable first
% solution point, our surgery inserts exactly NTSTINC mesh intervals such
% that the mesh of the excursion from the equilibrium remains unchanged
% (interpolation would introduce large errors), note also that the new
% period obtained here is significantly larger than above.

fac     = 50;
NTSTINC = 100;
NTST    = (numel(sol1.t)-1)/4; % NCOL = 4
PTINC   = NTSTINC*4+2; % NCOL = 4.
T       = sol1.t(end);
scale   = fac*NTSTINC;
ti = linspace(sol1.t(idx,1), T*(scale-1)+sol1.t(idx+1,1), PTINC);
xi = zeros(PTINC,3);
for i=1:3
  xi(:,i) = linspace(sol1.x(idx,i), sol1.x(idx+1,i), PTINC);
end
t0 = [sol1.t(1:idx-1,1) ; ti' ; T*(scale-1)+sol1.t(idx+2:end,1)];
x0 = [sol1.x(1:idx-1,:) ; xi ; sol1.x(idx+2:end,:)];
p0 = sol1.p;

prob = coco_prob();
prob = coco_set(prob, 'coll', 'NTST', NTST+NTSTINC, 'NTSTMX', 200);
% important detail: set TOL to something high while TOLINC and TOLDEC are
% set to values that force appropriate re-meshing --> MXCL does not occur
% while we obtain good mesh
prob = coco_set(prob, 'coll', 'TOL', 1.0e-1, 'TOLINC', 5.0e-5, 'TOLDEC', 1.0e-5);
prob = po_isol2orb(prob, '', @marsden, t0, x0, {'p1' 'p2'}, p0);

prob = coco_xchg_pars(prob, 'p2', 'po.period'); % Exchange parameters to fix period and free p2
prob = coco_set(prob, 'cont', 'NAdapt', 100); % Number of remesh-correct cycles
coco(prob, 'run2b', [], 0, {'p1' 'p2' 'po.period' 'po.seg.coll.NTST' 'po.seg.coll.err'}); % Iterate adaptation map

% for plot of dynamics of adaptation map use
%   bd = coco_bd_read('run2b');
%   plot(coco_bd_col(bd,'PT'), coco_bd_col(bd,'po.seg.coll.err'))
%   plot(coco_bd_col(bd,'PT'), coco_bd_col(bd,'po.seg.coll.NTST'))

coco_use_recipes_toolbox % remove the coll_v3, coll_v5, coll_v6, po_v1, and po_v3 toolboxes from the search path
