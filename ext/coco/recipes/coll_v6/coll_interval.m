function int = coll_interval(NCOL, dim)
%COLL_INTERVAL   Initialize data defining a collocation interval.
%
% Initialize properties that are independent of the discretization order (=
% number of mesh intervals) and the discretization parameters (= the
% warping coefficients). Structure does not change during continuation.
%
% Identical to coll_v4.
% 
% INT = COLL_INTERVAL(NCOL, DIM)
%
% INT  - Data structure.
% NCOL - Degree of interpolating polynomials.
% DIM  - State-space dimension.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_interval.m 2839 2015-03-05 17:09:01Z fschild $

int.NCOL = NCOL; % Degree of polynomial interpolants
int.dim  = dim; % Degree of polynomial interpolants

[int.tc int.wt] = coll_nodes(NCOL); % Degree of polynomial interpolants
int.tm = linspace(-1, 1, NCOL+1)';  % Basepoint temporal mesh
pmap   = coll_L(int.tm, int.tc);
dmap   = coll_Lp(int.tm, int.tc);
cmap   = coll_Lm(int.tm);
int.W  = kron(pmap, eye(dim)); % Interpolation matrix
int.Wp = kron(dmap, eye(dim)); % Interpolation matrix
int.Wm = kron(cmap, eye(dim)); % Interpolation matrix

end

function [nds wts] = coll_nodes(m)
%COLL_NODES   Compute collocation nodes and integration weights.
%
% Uses eigenvalues and eigenvectors of Jacobi matrix.
%
% Identical to coll_v1.
%
% [NDS WTS] = COLL_NODES(M)
%
% NDS - Collocation nodes.
% WTS - Quadrature weights.
% M   - Polynomial degree.

n = (1:m-1)';
g = n.*sqrt(1./(4.*n.^2-1));
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
% Identical to coll_v1.
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
% Identical to coll_v1.
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

function A = coll_Lm(ts)
%COLL_LM   Compute highest order coefficient of Lagrange polynomial.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% Identical to coll_v3.
%
% A = COLL_LM(TS)
% 
% A  - Array of coefficients.
% TS - Array of basepoints.

q = numel(ts);

sj = repmat(reshape(ts, [1 q]), [q 1]);
sk = repmat(reshape(ts, [q 1]), [1 q]);
t1 = sj-sk;
idx = abs(t1)<=eps;
t1(idx) = 1;

A = 1./t1;
A = prod(A, 2)';

end
