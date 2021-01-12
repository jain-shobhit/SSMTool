function [sol data] = msbvp_read_solution(oid, run, lab)
%MSBVP_READ_SOLUTION   Read 'msbvp' solution and toolbox data from disk.
%
% Extract data structure associated with 'msbvp' toolbox instance
% identifier and construct solution structure using trajectory segments
% from individual 'coll' toolbox instances.
%
% [SOL DATA] = MSBVP_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'msbvp');
data = coco_read_solution(tbid, run, lab);

sol = cell(1, data.nsegs);
for i=1:data.nsegs
  segoid = coco_get_id(tbid, sprintf('seg%d', i)); % 'coll' object instance identifier
  sol{i} = coll_read_solution(segoid, run, lab); % Trajectory segment
end

end
