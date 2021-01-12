function prob = doedel_sol2het(prob, run, lab)
%DOEDEL_SOL2HET   Construct instance of doedel continuation problem from stored data.
%
% PROB = DOEDEL_SOL2HET(PROB, RUN, LAB)
%
% PROB - continuation problem structure.
% RUN  - run identifier (string).
% LAB  - solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: doedel_sol2het.m 2839 2015-03-05 17:09:01Z fschild $

prob = coll_sol2seg(prob, 'doedel1', run, lab); % Construct first instance of 'coll' toolbox
prob = coll_sol2seg(prob, 'doedel2', run, lab); % Construct second instance of 'coll' toolbox
prob = alg_sol2eqn(prob, 'doedel3', run, lab);  % Construct first instance of 'alg' toolbox
prob = alg_sol2eqn(prob, 'doedel4', run, lab);  % Construct second instance of 'alg' toolbox

[data chart] = coco_read_solution('evs', run, lab);
vec  = chart.x(data.vec_idx); % Extract stored eigenvector
lam  = chart.x(data.lam_idx); % Extract stored eigenvalue

[data chart] = coco_read_solution('bcs', run, lab);
eps  = chart.x(data.eps_idx); % Extract stored deviations from equilibria
th   = chart.x(data.th_idx);  % Extract stored angle relative to horizontal

prob = doedel_close_het(prob, eps, th, vec, lam); % Apply gluing conditions

end
