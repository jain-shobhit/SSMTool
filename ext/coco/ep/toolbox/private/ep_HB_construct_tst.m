function [prob, data] = ep_HB_construct_tst(prob, data, sol)
%EP_HB_CONSTRUCT_TST   Add HB codim-2 test functions.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_HB_construct_tst.m 2869 2015-07-30 18:59:50Z hdankowicz $

data = init_data(data, sol);

ep = data.ep;

if any([ep.BTP])
  fid  = data.hb_tst.fid;
  pids = coco_get_id(fid, {'BTP'});
  
  if ep.BTP
    uidx = coco_get_func_data(prob, data.ep_hb.fid, 'uidx');
    prob = coco_add_pars(prob, fid, uidx(data.ep_hb.k_idx), pids, 'regular');
    prob = coco_add_event(prob, 'BTP', pids, 0);
  end
end

end

function data = init_data(data, ~)
%INIT_DATA   Initialize test function data.

tst.fid = coco_get_id(data.ep_eqn.fid, 'test');
data.hb_tst = tst;

end
