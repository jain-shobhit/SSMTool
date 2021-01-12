function data = coll_init_data(data, x0, p0)
%COLL_INIT_DATA   Initialize toolbox data for an instance of 'coll'.
%
% Populate remaining fields of the toolbox data structure used by 'coll'
% function objects.
%
% DATA = COLL_INIT_DATA(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for discretized trajectory.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_init_data.m 2839 2015-03-05 17:09:01Z fschild $

NTST = data.coll.NTST; % Number of mesh intervals
NCOL = data.coll.NCOL; % Degree of polynomial interpolants
dim  = size(x0, 2);    % State-space dimension
pdim = numel(p0);      % Number of problem parameters

data.dim  = dim;       % State-space dimension
data.pdim = pdim;      % Number of problem parameters

bpnum  = NCOL+1;            % Number of basepoints per interval
bpdim  = dim*(NCOL+1);      % Number of basepoint values per interval
xbpnum = (NCOL+1)*NTST;     % Number of basepoints
xbpdim = dim*(NCOL+1)*NTST; % Number of basepoint values
cndim  = dim*NCOL;          % Number of collocation node values per interval
xcnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions
cntnum = NTST-1;            % Number of internal boundaries
cntdim = dim*(NTST-1);      % Number of continuity conditions

data.xbp_idx = (1:xbpdim)'; % Index array for basepoint values
data.T_idx   = xbpdim+1;    % Index for interval length
data.p_idx   = xbpdim+1+(1:pdim)'; % Index array for problem parameters
data.tbp_idx = setdiff(1:xbpnum, 1+bpnum*(1:cntnum))'; % Index array without duplication of internal boundaries
data.x_shp   = [dim xcnnum]; % Shape for vectorization
data.xbp_shp = [dim xbpnum]; % Shape for vectorization
data.p_rep   = [1 xcnnum];   % Shape for vectorization

tm = linspace(-1, 1, bpnum)';
t  = repmat((0.5/NTST)*(tm+1), [1 NTST]);
t  = t+repmat((0:cntnum)/NTST, [bpnum 1]);
data.tbp    = t(:)/t(end);   % Temporal mesh with duplication

data.x0_idx = (1:dim)';            % Index array for trajectory end point at t=0
data.x1_idx = xbpdim-dim+(1:dim)'; % Index array for trajectory end point at t=1

[tc wts]    = coll_nodes(NCOL); % Collocation nodes and quadrature weights
wts         = repmat(wts, [dim NTST]);
data.wts1   = wts(1,:);                          
data.wts2   = spdiags(wts(:), 0, xcndim, xcndim);

pmap        = coll_L(tm, tc);
dmap        = coll_Lp(tm, tc);
rows        = reshape(1:xcndim, [cndim NTST]);
rows        = repmat(rows, [bpdim 1]);
cols        = repmat(1:xbpdim, [cndim 1]);
W           = repmat(kron(pmap, eye(dim)), [1 NTST]); % Interpolation matrix
Wp          = repmat(kron(dmap, eye(dim)), [1 NTST]); % Interpolation matrix
data.W      = sparse(rows, cols, W);
data.Wp     = sparse(rows, cols, Wp);

data.dxrows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);  % Index array for vectorization
data.dxcols = repmat(1:xcndim, [dim 1]);                         % Index array for vectorization
data.dprows = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]); % Index array for vectorization
data.dpcols = repmat(1:pdim, [dim xcnnum]);                      % Index array for vectorization

temp        = reshape(1:xbpdim, [bpdim NTST]);
Qrows       = [1:cntdim 1:cntdim];
Qcols       = [temp(1:dim, 2:end) temp(cndim+1:end, 1:end-1)];
Qvals       = [ones(cntdim,1) -ones(cntdim,1)];
data.Q      = sparse(Qrows, Qcols, Qvals, cntdim, xbpdim); % Jacobian of continuity conditions
data.dTpcnt = sparse(cntdim, 1+pdim);

end

function [nds wts] = coll_nodes(m)
%COLL_NODES   Compute collocation nodes and integration weights.
%
% Uses eigenvalues and eigenvectors of Jacobi matrix.
%
% [NDS WTS] = COLL_NODES(M)
%
% NDS - Collocation nodes.
% WTS - Quadrature weights.
% M   - Polynomial degree.

n = (1:m-1)';
g = n.*sqrt(1./(4*n.^2-1));
J = -diag(g,1)-diag(g,-1);

[w x] = eig(J);
nds   = diag(x);
wts   = 2*w(1,:).^2;

end

function A = coll_L(ts, tz)
%COLL_L   Evaluation of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_L(TS, TZ)
% 
% A  - Array of interpolated values.
% TS - Array of basepoints.
% TZ - Array of interpolation points.

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1]), [1 q q]);
sj = repmat(reshape(ts, [1 q 1]), [p 1 q]);
sk = repmat(reshape(ts, [1 1 q]), [p q 1]);

t1 = zi-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

A = prod(t1./t2, 3);

end

function A = coll_Lp(ts, tz)
%COLL_LP   Evaluation of derivative of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_LP(TS, TZ)
% 
% A  - Array of interpolated values
% TS - Array of basepoints
% TZ - Array of interpolation points

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1]), [1 q q q]);
sj = repmat(reshape(ts, [1 q 1 1]), [p 1 q q]);
sk = repmat(reshape(ts, [1 1 q 1]), [p q 1 q]);
sl = repmat(reshape(ts, [1 1 1 q]), [p q q 1]);

t3 = sj(:,:,:,1)-sk(:,:,:,1);
t4 = zi-sl;
t5 = sj-sl;

idx1 = find(abs(t5)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(sk-sl)<=eps);
t5(union(idx1, idx3)) = 1;
t4(union(idx1, idx3)) = 1;
t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

A = sum(t3.*prod(t4./t5, 4), 3);

end
