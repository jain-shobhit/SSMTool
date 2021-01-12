function [prob, data] = coll_construct_err(prob, data)
%COLL_CONSTRUCT_ERR   Add error test functions.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_err.m 3157 2019-12-16 23:36:58Z hdankowicz $

seg  = data.coll_seg;
fid  = seg.fid;

prob = coco_add_chart_data(prob, fid, struct(), struct());
uidx = coco_get_func_data(prob, fid, 'uidx');
efid = coco_get_id(fid, {'err' 'err_TF'});
prob = coco_add_func(prob, efid{1}, @coll_err, data, ...
  'regular', efid, 'uidx', uidx(seg.maps.xbp_idx), ...
  'remesh', @coll_err_remesh, 'passChart');
if data.coll.MXCL
  prob = coco_add_event(prob, 'MXCL', 'MX', efid{2}, '>', 1);
end
end

function [data, chart, y] = coll_err(prob, data, chart, u) %#ok<INUSL>

pr = data.pr;
seg = pr.coll_seg;
fid = seg.fid;

cdata = coco_get_chart_data(chart, fid);
if isfield(cdata, 'err')
  y = cdata.err;
else
  int  = seg.int;
  maps = seg.maps;
  
  cp = reshape(maps.Wm*u, [int.dim maps.NTST]);
  y  = maps.wn*max(sqrt(sum(cp.^2,1)));
  y  = [y; y/pr.coll.TOL];
  cdata.err = y;
  chart = coco_set_chart_data(chart, fid, cdata);
end

end

function [prob, stat, xtr] = coll_err_remesh(prob, data, chart, ub, Vb) %#ok<INUSD>

seg  = data.coll_seg;
maps = seg.maps;

xtr  = [];
uidx = coco_get_func_data(prob, seg.fid, 'uidx');
prob = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
stat = 'success';

end
