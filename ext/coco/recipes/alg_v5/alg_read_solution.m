function [sol data] = alg_read_solution(oid, run, lab)
%ALG_READ_SOLUTION   Read 'alg' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'alg' toolbox instance
% identifier from solution file and construct solution structure.
%
% Differs from alg_v4 by providing support for unique toolbox instance
% identifiers.
%
% [SOL DATA] = ALG_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - object instance identifier (string).
% RUN  - run identifier (string).
% LAB  - solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: alg_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid         = coco_get_id(oid, 'alg'); % Create toolbox instance identifier
[data chart] = coco_read_solution(tbid, run, lab);
sol.x = chart.x(data.x_idx);
sol.p = chart.x(data.p_idx);
sol.u = [sol.x; sol.p];

end
