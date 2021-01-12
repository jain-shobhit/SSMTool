function [prob, status, xtr, ftr] = adjt_int_remesh(prob, data, chart, lb, Vlb)  %#ok<INUSL>

cid = coco_get_id(data.oid, 'coll');
[fdata, axidx] = coco_get_adjt_data(prob, cid, 'data', 'axidx');
data = adjt_int_init_data(fdata, data.oid);
opt  = data.coll_opt;
aidx = axidx([opt.xcn_idx; opt.T_idx]);

xtr    = [];
ftr    = 1;
prob   = coco_change_adjt(prob, data, 'aidx', aidx, 'l0', lb, 'vecs', Vlb);
status = 'success';

end
