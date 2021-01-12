function [prob, status, xtr] = int_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

cid = coco_get_id(data.oid, 'coll');
[fdata, uidx] = coco_get_func_data(prob, cid, 'data', 'uidx');
data = int_init_data(fdata, data.oid);
maps = data.coll_seg.maps;
uidx = uidx([maps.xbp_idx; maps.T_idx]);

xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end
