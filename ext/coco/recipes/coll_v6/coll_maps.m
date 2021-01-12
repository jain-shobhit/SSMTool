function maps = coll_maps(int, NTST, pdim)
%COLL_MAPS   Initialise data depending on order but not on mesh distribution.
%
% Initialize properties that are independent of the discretization order (=
% number of mesh intervals) and the discretization parameters (= the
% mesh distribution).
%
% Differs from coll_v5 by providing support for varying discretization
% order.
%
% MAPS = COLL_MAPS(INT, NTST, PDIM)
%
% MAPS - Data structure.
% INT  - Data structure.
% NTST - Number of mesh intervals.
% PDIM - Number of problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_maps.m 2839 2015-03-05 17:09:01Z fschild $

NCOL = int.NCOL; % Degree of polynomial interpolants
dim  = int.dim;  % State-space dimension

maps.NTST = NTST; % Number of mesh intervals
maps.pdim = pdim; % Number of problem parameters

bpnum  = NCOL+1;            % Number of basepoints per interval
bpdim  = dim*(NCOL+1);      % Number of basepoint values per interval
xbpnum = (NCOL+1)*NTST;     % Number of basepoints
xbpdim = dim*(NCOL+1)*NTST; % Number of basepoint values
cndim  = dim*NCOL;          % Number of collocation node values per interval
xcnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions
cntnum = NTST-1;            % Number of internal boundaries
cntdim = dim*(NTST-1);      % Number of internal boundary values

maps.xbp_idx = (1:xbpdim)'; % Index array for basepoint values
maps.T_idx   = xbpdim+1;    % Index for interval length
maps.p_idx   = xbpdim+1+(1:pdim)';      % Index array for problem parameters
maps.Tp_idx  = [maps.T_idx; maps.p_idx];
maps.tbp_idx = setdiff(1:xbpnum, 1+bpnum*(1:cntnum))';      % Index array for problem parameters
maps.x_shp   = [dim xcnnum]; % Shape for vectorization
maps.xbp_shp = [dim xbpnum]; % Shape for vectorization
maps.p_rep   = [1 xcnnum];   % Shape for vectorization

maps.xtr     = [maps.xbp_idx; maps.Tp_idx]; % Shape for vectorization
maps.xtr(dim+1:end-dim-pdim-1) = 0; % Boundary elements, interval length, and parameters
maps.xtrend  = maps.xtr(end-dim-pdim:end); % Right boundary, interval length, and parameters

rows         = reshape(1:xcndim, [cndim NTST]);
rows         = repmat(rows, [bpdim 1]);
cols         = repmat(1:xbpdim, [cndim 1]);
W            = repmat(int.W, [1 NTST]);  % Interpolation matrix
Wp           = repmat(int.Wp, [1 NTST]); % Interpolation matrix
maps.W       = sparse(rows, cols, W);
maps.Wp      = sparse(rows, cols, Wp);

temp         = reshape(1:xbpdim, [bpdim NTST]);
Qrows        = [1:cntdim 1:cntdim];
Qcols        = [temp(1:dim, 2:end) temp(cndim+1:end, 1:end-1)];
Qvals        = [ones(cntdim,1) -ones(cntdim,1)];
maps.Q       = sparse(Qrows, Qcols, Qvals, cntdim, xbpdim); % Jacobian of continuity conditions
maps.Qnum    = cntdim; % Number of rows for Jacobian of continuity conditions with respect to T and p.

maps.dxrows  = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);  % Index array for vectorization
maps.dxcols  = repmat(1:xcndim, [dim 1]);                         % Index array for vectorization
maps.dprows  = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]); % Index array for vectorization
maps.dpcols  = repmat(1:pdim, [dim xcnnum]);                      % Index array for vectorization

maps.x0_idx  = (1:dim)';            % Index array for trajectory end point at t=0
maps.x1_idx  = xbpdim-dim+(1:dim)'; % Index array for trajectory end point at t=1

rows         = reshape(1:dim*NTST, [dim NTST]);
rows         = repmat(rows, [bpdim 1]);
cols         = repmat(1:xbpdim, [dim 1]);
Wm           = repmat(int.Wm, [1 NTST]); % Index array for trajectory end point at t=1
maps.Wm      = sparse(rows, cols, Wm);
x            = linspace(int.tm(1), int.tm(2), 51);
y            = arrayfun(@(x) prod(x-int.tm), x);
maps.wn      = max(abs(y)); % Bound on mesh products

end
