function data = adjt_int_init_data(fdata, oid)

data.coll_seg  = fdata.coll_seg;
data.coll_opt  = fdata.coll_opt;
data.xbp_idx   = data.coll_seg.maps.xbp_idx;
data.T_idx     = data.xbp_idx(end) + 1;
data.ghan      = @ghan;
data.ghan_dx   = @ghan_dx;
data.ghan_dxdx = @ghan_dxdx;

seg  = fdata.coll_seg;
maps = seg.maps;
int  = seg.int;

NCOL = int.NCOL;
NTST = maps.NTST;
xdim = int.dim;

cndim  = NCOL*xdim;
bpdim  = xdim*(NCOL+1);
xbpdim = NTST*(NCOL+1)*xdim;
xcnnum = NTST*NCOL;
xcndim = NTST*NCOL*xdim;
addim  = xcndim+1;

% Derivative of (T/2N)*gxcn with respect to xbp:
rows = 1 + xcnnum*repmat(0:xdim-1, [xdim 1]);
rows = repmat(rows(:), [1 xcnnum]) + repmat(0:xcnnum-1, [xdim^2 1]);
opt.gdxdxrows1 = rows;
cols = 1:xcndim;
opt.gdxdxcols1 = repmat(reshape(cols, [xdim xcnnum]), [xdim 1]);

opt.gdxdxrows2 = ones(cndim*xbpdim, 1);
cols = 1 + xdim*(0:NCOL-1);
cols = repmat(cols(:), [1 xdim]) + repmat(0:xdim-1, [NCOL 1]);
cols = repmat(cols(:), [1 bpdim]) + addim*repmat(0:bpdim-1, [cndim 1]);
cols = repmat(cols(:), [1 NTST]) + ...
  (cndim+addim*bpdim)*repmat(0:NTST-1, [cndim*bpdim 1]);
opt.gdxdxcols2 = cols;

idx = 1:NCOL;
idx = repmat(idx(:), [1 xdim]) + xcnnum*repmat(0:xdim-1, [NCOL 1]);
idx = repmat(idx(:), [1 bpdim]) + xcndim*repmat(0:bpdim-1, [cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  (NCOL+cndim*xbpdim)*repmat(0:NTST-1, [cndim*bpdim 1]);
opt.gdxdxidx = idx;

% Derivative of (T/2N)*gxcn with respect to T:
opt.gdxdTrows = ones(xcndim,1);
opt.gdxdTcols = addim*xbpdim  + (1:xcndim)';

% Derivative of (1/2N)*w*g' with respect to xbp:
opt.gdTdxrows = ones(xbpdim,1);
opt.gdTdxcols = xcndim + 1 + addim*(0:xbpdim-1)';

opt.dJrows = 1;
opt.dJcols = addim*(xbpdim+1);

data.int_opt = opt;
data.oid = oid;

data = coco_func_data(data);

end
