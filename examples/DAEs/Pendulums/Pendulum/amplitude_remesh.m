function [prob, status, xtr] = amplitude_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[colldata, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
maps = colldata.coll_seg.maps;
xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
status = 'success';

end