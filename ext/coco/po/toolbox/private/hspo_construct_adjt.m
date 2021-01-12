function prob = hspo_construct_adjt(prob, data, sol)
%HSPO_CONSTRUCT_ADJT   Add HSPO adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: po_construct_adjt.m 2872 2015-08-06 20:12:06Z hdankowicz $

data = init_data(data);
opt  = data.hspo_opt;
fid  = opt.fid;

[fdata, axidx]= coco_get_adjt_data(prob, data.bvid, 'data', 'axidx');
pfid = coco_get_id(fid, 'period');
dpar = coco_get_id('d', pfid);
aidx = axidx(fdata.bvp_bc.T_idx);
prob = coco_add_adjt(prob, pfid, dpar, 'aidx', aidx, ...
  'l0', sol.period_l0, 'tl0', sol.period_tl0);

end

function data = init_data(data)
%INIT_DATA   Initialize HSPO data with phase condition

opt.fid = coco_get_id(data.oid, 'hspo');
data.hspo_opt = opt;

end
