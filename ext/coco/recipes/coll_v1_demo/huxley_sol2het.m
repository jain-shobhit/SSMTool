function prob = huxley_sol2het(prob, run, lab)
%HUXLEY_SOL2HET   Construct instance of huxley continuation problem from stored data.
%
% PROB = HUXLEY_SOL2HET(PROB, RUN, LAB)
%
% PROB - Continuation problem structure.
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: huxley_sol2het.m 2839 2015-03-05 17:09:01Z fschild $

prob = coll_sol2seg(prob, 'huxley1', run, lab); % Construct first instance of 'coll' toolbox
prob = coll_sol2seg(prob, 'huxley2', run, lab); % Construct first instance of 'coll' toolbox

[data chart] = coco_read_solution('bcs', run, lab);
epsv = chart.x(data.eps_idx); % Extract stored deviations of heteroclinic trajectory end points from corresponding equilibria

prob = huxley_close_het(prob, epsv); % Apply gluing conditions

end
