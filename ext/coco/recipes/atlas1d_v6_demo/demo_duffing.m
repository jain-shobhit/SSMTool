coco_use_recipes_toolbox atlas1d_v6 po_v1 coll_v1 % Add atlas1d_v6 atlas algorithm and po_v1 and coll_v1 toolboxes to search path

%% Example 16.2

% The continuation problem structure encoded below corresponds to a family
% of 601 zero functions (480 collocation conditions, 116 continuity
% conditions, four periodic boundary conditions, one integral phase
% condition) in terms of 606 continuation variables (600 basepoint values,
% one interval length, and five problem parameters), a family of five
% monitor function that evaluate to the problem parameters, and the
% corresponding inactive continuation parameters 'la', 'al', 'eps', 'A',
% and 'om'. Its dimensional deficit equals 0. A frequency-response curve is
% obtained by releasing 'om' and allowing it to vary during continuation.

p0 = [0.2; 1; 1; 2.5; 1];
x0 = [0; 0; 1; 0];
[t0 x0] = ode45(@(t,x) duff(x,p0), [0 30], x0); % Let transients die out
[t0 x0] = ode45(@(t,x) duff(x,p0), [0 2*pi], x0(end,:)'); % An approximate periodic orbit
prob = coco_set(coco_prob(), 'coll', 'NTST', 30);
prob = po_isol2orb(prob, '', @duff, t0, x0, ...
  {'la' 'al' 'eps' 'A' 'om'}, p0);
prob = coco_set(prob, 'cont', 'atlas', @atlas_1d_min.create);
prob = coco_set(prob, 'cont', 'PtMX', 500, 'h', 0.5, 'almax', 30);
coco(prob, 'duffing', [], 1, 'om', [0.5 3.5]);

coco_use_recipes_toolbox % Remove atlas1d_v6 atlas algorithm and po_v1 and coll_v1 toolboxes from search path
