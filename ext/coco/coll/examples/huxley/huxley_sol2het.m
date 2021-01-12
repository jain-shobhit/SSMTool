function prob = huxley_sol2het(prob, run, lab)
%HUXLEY_SOL2HET   Construct instance of huxley continuation problem from stored data.
%
% PROB = HUXLEY_SOL2HET(PROB, RUN, LAB)
%
% PROB - Continuation problem structure.
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: huxley_sol2het.m 2862 2015-07-26 22:11:32Z hdankowicz $

prob = ode_coll2coll(prob, 'huxley1', run, lab); % Construct first instance of 'coll' toolbox
prob = ode_coll2coll(prob, 'huxley2', run, lab); % Construct first instance of 'coll' toolbox

[data, chart] = coco_read_solution('bcs', run, lab);
devs = chart.x(data.dev_idx); % Extract stored deviations of heteroclinic trajectory end points from corresponding equilibria

prob = huxley_close_het(prob, devs); % Apply gluing conditions

end
