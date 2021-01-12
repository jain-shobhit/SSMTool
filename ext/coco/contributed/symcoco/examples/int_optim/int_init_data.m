%% Initialize |data| for objective functional of integral type
% This function uses symbolic-toolbox generated derivatives of the
% integrand. See <demo.html> and <gen_sym_int_optim.html> for other parts
% of the demo, and |coco_folder/po/examples/int_optim| and PO-Tutorial.pdf
% for details about the demo.
%%
function data = int_init_data(fdata, oid)
data.coll_seg = fdata.coll_seg;
data.xbp_idx  = data.coll_seg.maps.xbp_idx;
data.T_idx    = data.xbp_idx(end) + 1;
%% Shortcut G for defining functional integrand and its derivatives
G=sco_gen(@sym_g);
data.ghan     = G('');
data.ghan_dx  = G('x');
data.oid      = oid;

data = coco_func_data(data);

end
