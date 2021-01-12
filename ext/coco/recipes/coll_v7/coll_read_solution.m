function [sol data] = coll_read_solution(oid, run, lab)
%COLL_READ_SOLUTION   Read 'coll' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'coll' toolbox instance
% identifier from solution file and construct solution structure.
%
% Identical to coll_v4.
%
% [SOL DATA] = COLL_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_read_solution.m 2839 2015-03-05 17:09:01Z fschild $

tbid         = coco_get_id(oid, 'coll');
[data chart] = coco_read_solution(tbid, run, lab);
maps         = data.maps;
mesh         = data.mesh;

sol.t = mesh.tbp(maps.tbp_idx)*chart.x(maps.T_idx);
xbp   = reshape(chart.x(maps.xbp_idx), maps.xbp_shp)';
sol.x = xbp(maps.tbp_idx,:);
sol.p = chart.x(maps.p_idx);

end
