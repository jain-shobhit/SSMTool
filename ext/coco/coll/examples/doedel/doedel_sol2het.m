function prob = doedel_sol2het(prob, run, lab)
%DOEDEL_SOL2HET   Construct instance of doedel continuation problem from stored data
%
% PROB = DOEDEL_SOL2HET(PROB, RUN, LAB)
%
% PROB - continuation problem structure.
% RUN  - run identifier (string).
% LAB  - solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: doedel_sol2het.m 2896 2015-10-05 20:08:38Z hdankowicz $

% Reconstruct the two ‘coll’ instances.
prob = ode_coll2coll(prob, 'doedel1', run, lab);
prob = ode_coll2coll(prob, 'doedel2', run, lab);

% Reconstruct the two ‘ep’ instances. A single-vector variational
% problem is automatically included with the second instance.
prob = ode_ep2ep(prob, 'doedel3', run, lab);  % Construct first instance of 'ep' toolbox
prob = ode_ep2ep(prob, 'doedel4', run, lab);  % Construct second instance of 'ep' toolbox

% Extract stored eigenvalue
[data, chart] = coco_read_solution('evs', run, lab);
lam  = chart.x(data.lam_idx);

 % Extract stored deviations from equilibria and angle relative to horizontal
[data, chart] = coco_read_solution('bcs', run, lab);
dev  = chart.x(data.dev_idx);
th   = chart.x(data.th_idx);

% Apply gluing conditions
prob = doedel_close_het(prob, dev, th, lam);

end
