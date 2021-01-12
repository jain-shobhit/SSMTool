function [sol data] = alg_read_solution(run, lab)
%ALG_READ_SOLUTION   Read 'alg' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'alg' toolbox instance
% identifier from solution file and construct solution structure.
%
% Differs from alg_v3 by returning data and solution structure instead of
% just arrays of problem variables and problem parameters. Implementation
% supports using alg_read_solution when constructing an instance of 'alg'
% from stored data.
% 
% [SOL DATA] = ALG_READ_SOLUTION(RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

[data chart] = coco_read_solution('alg', run, lab);
sol.x = chart.x(data.x_idx);
sol.p = chart.x(data.p_idx);
sol.u = [sol.x; sol.p];

end
