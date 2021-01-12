function [sol data] = po_read_solution(oid, run, lab)
%PO_READ_SOLUTION   Read 'po' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'coll' and 'po' toolbox
% instance identifier from solution file and construct solution structure.
%
% Identical to po_v1.
%
% [SOL DATA] = PO_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'po');
data = coco_read_solution(tbid, run, lab);

segoid = coco_get_id(tbid, 'seg');
sol    = coll_read_solution(segoid, run, lab); % Trajectory segment

end
