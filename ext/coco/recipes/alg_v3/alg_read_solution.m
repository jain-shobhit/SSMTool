function [x p] = alg_read_solution(run, lab)
%ALG_READ_SOLUTION   Read 'alg' solution from disk.
%
% Load data and chart structures associated with 'alg' toolbox instance
% identifier from solution file and extract arrays of problem variables and
% problem parameters.
% 
% [X P] = ALG_READ_SOLUTION(RUN, LAB)
%
% X   - Problem variables (array).
% P   - Problem parameters (array).
% RUN - Run identifier (string).
% LAB - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

[data chart] = coco_read_solution('alg', run, lab);
x = chart.x(data.x_idx);
p = chart.x(data.p_idx);

end
