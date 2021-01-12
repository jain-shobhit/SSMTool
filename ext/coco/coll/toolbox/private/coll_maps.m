function maps = coll_maps(int, NTST, pdim)

NCOL = int.NCOL;
dim  = int.dim;

maps.NTST = NTST;
maps.pdim = pdim;

bpnum  = NCOL+1;
bpdim  = dim*(NCOL+1);
xbpnum = (NCOL+1)*NTST;
xbpdim = dim*(NCOL+1)*NTST;
cndim  = dim*NCOL;
xcnnum = NCOL*NTST;
xcndim = dim*NCOL*NTST;
cntnum = NTST-1;
cntdim = dim*(NTST-1);

maps.fdim    = xcndim + cntdim;
maps.xbp_idx = (1:xbpdim)';
maps.T0_idx  = xbpdim+1;
maps.T_idx   = xbpdim+2;
maps.p_idx   = maps.T_idx+(1:pdim)';
maps.Tp_idx  = [maps.T0_idx; maps.T_idx; maps.p_idx];
maps.tbp_idx = setdiff(1:xbpnum, 1+bpnum*(1:cntnum))';
maps.x_shp   = [dim xcnnum];
maps.xbp_shp = [dim xbpnum];
maps.p_rep   = [1 xcnnum];

maps.xtr     = [maps.xbp_idx; maps.Tp_idx];
maps.xtr(dim+1:end-dim-pdim-2) = 0;
maps.xtrend  = maps.xtr(end-dim-pdim-1:end); % Right boundary, start time, interval length, and parameters

rows         = reshape(1:xcndim, [cndim NTST]);
rows         = repmat(rows, [bpdim 1]);
cols         = repmat(1:xbpdim, [cndim 1]);
W            = repmat(int.W, [1 NTST]);
Wp           = repmat(int.Wp, [1 NTST]);
maps.W       = sparse(rows, cols, W);
maps.Wp      = sparse(rows, cols, Wp);

temp         = reshape(1:xbpdim, [bpdim NTST]);
Qrows        = [1:cntdim 1:cntdim];
Qcols        = [temp(1:dim, 2:end) temp(cndim+1:end, 1:end-1)];
Qvals        = [ones(cntdim,1) -ones(cntdim,1)];
maps.Q       = sparse(Qrows, Qcols, Qvals, cntdim, xbpdim);
maps.Qnum    = cntdim;

maps.fdxrows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);
maps.fdxcols = repmat(1:xcndim, [dim 1]);
maps.fdprows = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]);
maps.fdpcols = repmat(1:pdim, [dim xcnnum]);
maps.gdxrows = repmat(1:xcnnum, [dim 1]);
maps.gdxcols = 1:xcndim;
maps.gdprows = repmat(1:xcnnum, [pdim 1]);
maps.gdpcols = repmat(1:pdim, [1, xcnnum]);

maps.x0_idx  = (1:dim)';
maps.x1_idx  = xbpdim-dim+(1:dim)';

rows         = reshape(1:dim*NTST, [dim NTST]);
rows         = repmat(rows, [bpdim 1]);
cols         = repmat(1:xbpdim, [dim 1]);
Wm           = repmat(int.Wm, [1 NTST]);
maps.Wm      = sparse(rows, cols, Wm);
x            = linspace(int.tm(1), int.tm(2), 51);
y            = arrayfun(@(x) prod(x-int.tm), x);
maps.wn      = max(abs(y));

end
