function data = var_coll_init_data(prob, data)
%VAR_COLL_INIT_DATA   Initialize toolbox data for an instance of 'varcoll'.
%
% Populate remaining fields of the toolbox data structure used by 'varcoll'
% function objects.
%
% DATA = VAR_COLL_INIT_DATA(PROB, DATA)
%
% DATA - Toolbox data structure.
% PROB - Continuation problem structure.

fdata = coco_get_func_data(prob, data.coll_id, 'data'); % Extract 'coll' instance data

NTST   = fdata.coll.NTST; % Number of mesh intervals
NCOL   = fdata.coll.NCOL; % Degree of interpolating polynomials
dim    = fdata.dim;       % State-space dimension
pdim   = fdata.pdim;      % Number of problem parameters
xbpnum = (NCOL+1)*NTST;   % Number of base points
xbpdim = dim*(NCOL+1)*NTST;   % Number of basepoint values for trajectory
ubpdim = dim^2*(NCOL+1)*NTST; % Number of basepoint values for fundamental solution
xcnnum = NCOL*NTST;       % Number of collocation nodes
xcndim = dim*NCOL*NTST;   % Number of collocation values for trajectory
ucndim = dim^2*NCOL*NTST; % Number of collocation values for fundamental solution

data.dim     = dim; % State-space dimension
data.M1_idx  = xbpdim-dim+(1:dim)'; % Index array for last point
data.ubp_idx = xbpdim+1+pdim+(1:ubpdim)'; % Index array for basepoint values for fundamental solution
data.u_shp   = [xbpdim dim]; % Shape for vectorization
data.R       = [speye(dim) sparse(dim, xbpdim-dim)];
data.Id      = eye(dim);
data.jac     = [sparse(dim^2*NTST, xbpdim+1+pdim) ...
               [kron(eye(dim), fdata.Q); kron(eye(dim), data.R)]];

rows = reshape(1:ucndim, [dim^2 xcnnum]);
data.dxdxrows1 = repmat(rows, [dim 1]);
data.dxdxcols1 = repmat(1:xcndim, [dim^2 1]);

rows = reshape(1:xbpdim*xcndim, [dim xbpnum*xcndim]);
data.dxdxrows2 = repmat(rows, [dim 1]);
data.dxdxcols2 = repmat(1:xcndim, [dim xbpdim]);

rows = reshape(1:ucndim, [xcndim dim]);
data.dxdxrows3 = repmat(rows, [xbpdim 1]);
data.dxdxcols3 = repmat(1:xbpdim, [xcndim dim]);

rows = reshape(1:xcndim, [dim xcnnum]);
data.dxdprows  = repmat(rows, [dim*pdim 1]);
cols = permute(reshape(1:xcndim*pdim, [dim xcnnum pdim]), [1 3 2]);
data.dxdpcols  = repmat(cols(:)', [dim 1]);
data.dxdp_shp  = [ucndim pdim];

end
