function int = coll_interval(NCOL, dim)

int.NCOL = NCOL;
int.dim  = dim;

[int.tc, int.wt] = coll_nodes(NCOL);
int.tm = linspace(-1, 1, NCOL+1)';
pmap   = coll_L(int.tm, int.tc);
dmap   = coll_Lp(int.tm, int.tc);
mmap   = coll_Lm(int.tm);
int.W  = kron(pmap, eye(dim));
int.Wp = kron(dmap, eye(dim));
int.Wm = kron(mmap, eye(dim));

end

function [nds, wts] = coll_nodes(m)
% implementation of Golub's method

n = (1:m-1)';
g = n.*sqrt(1./(4.*n.^2-1));
J = -diag(g,1)-diag(g,-1);

[w, x] = eig(J);
nds   = diag(x);
wts   = 2*w(1,:).^2;

end

function A = coll_L(ts, tz)

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

q = numel(ts);

sj = repmat(reshape(ts, [1 q]), [q 1]);
sk = repmat(reshape(ts, [q 1]), [1 q]);
t1 = sj-sk;
idx = abs(t1)<=eps;
t1(idx) = 1;

A = 1./t1;
A = prod(A, 2)';

end
