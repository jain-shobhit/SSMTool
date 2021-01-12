function [data sol] = calcvar_system(data, t0, x0, p0)
%CALCVAR_SYSTEM   Initialise function data and solution guess for quadrature discretization of catenary problem.
%
% [DATA SOL] = CALCVAR_SYSTEM(DATA, T0, X0, P0)
%
% DATA - Function data structure.
% SOL  - Initial solution guess.
% T0   - Array of temporal mesh points.
% X0   - Array of state vectors at mesh points.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: calcvar_system.m 2839 2015-03-05 17:09:01Z fschild $

NTST  = data.NTST; % Number of mesh intervals
NCOL  = data.NCOL; % Degree of interpolating polynomials

data.tk = linspace(-1, 1, NCOL+1); % Basepoint mesh on each interval
[data.th wt] = coll_nodes(NCOL);   % Collocation nodes and quadrature weights

tt = (0:NTST-1)' / NTST;
t  = repmat(data.tk, [NTST 1]);
t  = (t+1) * 0.5/NTST;
t  = t + repmat(tt, [1 NCOL+1]);
t  = reshape(t', [NTST*(NCOL+1) 1]);
t  = t./t(end);
data.tbp = t; % Basepoint mesh with duplication

T0   = t0(end) - t0(1);
t0   = (t0 - t0(1)) / T0; % Rescaling of independent variable
x0   = interp1(t0, x0, t)'; % Basepoint values
x0   = reshape(x0, [(NCOL+1)*NTST 1]);
pdim = numel(p0);
sol.u = [x0; p0];

wt   = repmat(reshape(wt', [NCOL 1]), [NTST 1]); % Quadrature weights
data.wt = spdiags(wt,0,NTST*NCOL,NTST*NCOL);

data.x_idx  = 1:(NCOL+1)*NTST;          % Index array for basepoint values
data.p_idx  = (NCOL+1)*NTST + (1:pdim); % Index array for problem parameter

intidx         = reshape(1:NTST*(NCOL+1), [NCOL+1, NTST]);
data.fint1_idx = intidx(end,1:end-1); % Index array of right end points for intervals 1..N-1
data.fint2_idx = intidx(1,2:end);     % Index array for left end points for intervals 2..N
intidx         = intidx(2:end,:);     % Skip left end points
data.fint3_idx = intidx(1:end-1);     % Index array for interior points without duplication

temp        = reshape(1:(NCOL+1)*NTST, [NCOL+1 NTST]);
ipidx       = reshape(temp(1, 2:end),[1 NTST-1]);
epidx       = reshape(temp(NCOL+1:end, 1:end-1),[1 NTST-1]);

rows        = [1:(NTST-1) 1:(NTST-1)];
cols        = [ipidx epidx];
vals        = [ones(1,NTST-1) -ones(1,NTST-1)];
data.Q      = sparse(rows, cols, vals, NTST-1, (NCOL+1)*NTST); % Coefficient matrix for continuity conditions

rows        = reshape(1:NCOL*NTST, [NCOL NTST]);
rows        = repmat(rows, [NCOL+1 1]);
rows        = reshape(rows, [(NCOL+1)*NCOL*NTST 1]);
cols        = repmat(1:(NCOL+1)*NTST, [NCOL 1]);
cols        = reshape(cols, [(NCOL+1)*NCOL*NTST 1]);

pmap        = coll_L(data.tk, data.th);
dmap        = coll_Lp(data.tk, data.th);
W           = repmat(pmap, [1 NTST]);
W           = reshape(W, [(NCOL+1)*NCOL*NTST 1]);
Wp          = repmat(dmap, [1 NTST]);
Wp          = reshape(Wp, [(NCOL+1)*NCOL*NTST 1]);
data.W      = sparse(rows, cols, W);  % Interpolation matrix
data.Wp     = sparse(rows, cols, Wp); % Interpolation matrix

end

function [x w] = coll_nodes(n)
%COLL_NODES   Compute collocation nodes and integration weights.
%
% Uses eigenvalues and eigenvectors of Jacobi matrix.
%
% [NDS WTS] = COLL_NODES(M)
%
% NDS - Collocation nodes.
% WTS - Quadrature weights.
% M   - Polynomial degree.

nn = 1:n-1;
ga = -nn.*sqrt(1./(4.*nn.^2-1));
J  = zeros(n,n);
J(sub2ind([n n], nn, nn+1)) = ga;
J(sub2ind([n n], nn+1, nn)) = ga;

[w,x] = eig(J);

x = diag(x);
w = 2*w(1,:)'.^2;

end

function A = coll_L(tk, th)
%COLL_L   Evaluation of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_L(TS, TZ)
% 
% A  - Array of interpolated values.
% TS - Array of basepoints.
% TZ - Array of interpolation points.

q = length(tk);
p = length(th);

ti = reshape(tk, [1 1 q]);
tk = reshape(tk, [1 q 1]);
th = reshape(th, [p 1 1]);

ti = repmat(ti, [p q 1]);
tk = repmat(tk, [p 1 q]);
th = repmat(th, [1 q q]);

tki = tk-ti;
thi = th-ti;

idx = find(abs(tki)<=eps);

thi(idx) = 1;
tki(idx) = 1;

A = thi./tki;
A = prod(A, 3);

A = reshape(A, [p q]);

end

function A = coll_Lp(tk, th)
%COLL_LP   Evaluation of derivative of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_LP(TS, TZ)
% 
% A  - Array of interpolated values
% TS - Array of basepoints
% TZ - Array of interpolation points

q = length(tk);
p = length(th);

tj = reshape(tk, [1 1 q]);
tj = repmat (tj, [p q 1]);

ti = reshape(tk, [1 1 1 q]);
tk = reshape(tk, [1 q 1 1]);
th = reshape(th, [p 1 1 1]);

ti = repmat(ti, [p q q 1]);
tk = repmat(tk, [p 1 q q]);
th = repmat(th, [1 q q q]);

tki = tk-ti;
tkj = tk(:,:,:,1)-tj;
thi = th-ti;

idx1 = find(abs(tki)<=eps);
idx2 = find(abs(tkj)<=eps);
idx3 = find(abs(repmat(tj, [1 1 1 q])-ti)<=eps);

tki(idx1) = 1;
tki(idx3) = 1;
thi(idx1) = 1;
thi(idx3) = 1;

tkj(idx2) = 1;
tkj       = 1.0 ./ tkj;
tkj(idx2) = 0;

A = thi ./ tki;
A = prod(A, 4);

A = reshape(A, [p q q]);
A = tkj .* A;
A = sum(A, 3);

A = reshape(A, [p q]);

end
