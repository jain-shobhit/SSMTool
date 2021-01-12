function [sol data] = compalg_read_solution(oid, run, lab)
%COMPALG_READ_SOLUTION   Read 'compalg' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'compalg' toolbox
% instance identifier from solution file and construct solution structure.
%
% Differs from compalg_v1 by support for using compalg_read_solution when
% constructing an instance of 'compalg' from stored data.
%
% [SOL DATA] = COMPALG_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID - object instance identifier (string).
% RUN - run identifier (string).
% LAB - solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: compalg_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'compalg'); % Create toolbox instance identifier
data = coco_read_solution(tbid, run, lab);

sol(data.neqs) = struct('x', [], 'p', []); % Allocate array with fields
for i=1:data.neqs
  soid = coco_get_id(tbid, sprintf('eqn%d', i)); % Create 'alg' object instance identifier
  algsol = alg_read_solution(soid, run, lab);
  sol(i).x = algsol.x;
  sol(i).p = algsol.p;
end

end
