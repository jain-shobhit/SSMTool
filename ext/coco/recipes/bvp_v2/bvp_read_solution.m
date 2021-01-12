function [sol data] = bvp_read_solution(oid, run, lab)
%BVP_READ_SOLUTION   Read 'bvp' solution and toolbox data from disk.
%
% Extract data structure associated with 'bvp' toolbox instance
% identifier and construct solution structure using trajectory segment
% from 'coll' toolbox instance.
%
% Identical to bvp_v1.
%
% [SOL DATA] = BVP_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: bvp_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'bvp');
data = coco_read_solution(tbid, run, lab);

segoid = coco_get_id(tbid, 'seg');
sol    = coll_read_solution(segoid, run, lab); % Trajectory segment

end
