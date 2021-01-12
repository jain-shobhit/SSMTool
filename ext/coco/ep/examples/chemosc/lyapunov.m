function [data, y] = lyapunov(prob, data, u) %#ok<INUSL>
%LYAPUNOV Monitor the first Lyapunov coefficient
%
% Computation of first Lyapunov coefficient for vector field quadratic in
% the state. The sign of the first Lyapunov coefficient determines the sub-
% or supercritical nature of a Hopf bifurcation.

x = u(data.x_idx);
p = u(data.p_idx);

% compute eigenvector
A = bykov_dx(x,p);
[X, D]   = eig(A);
[m, idx] = min(abs(real(diag(D))));
v  = X(:,idx);
om = imag(D(idx,idx));
vb = conj(v);
% if neutral saddle
if m>1e-6
  y = NaN;
  return
end

% compute adjoint eigenvector
[X, D]   = eig(A');
[m, idx] = min(abs(real(diag(D)))); %#ok<ASGLU>
w = X(:,idx);

if om*imag(D(idx,idx))>0
  w = conj(w);
end
w = w/conj(w'*v);

% compute tensor products
B  = bykov_dxdx(x,p);
B1 = zeros(numel(x),1);
B3 = zeros(numel(x),1);
for i=1:numel(x)
  Bmat  = reshape(B(i,:,:),[numel(x),numel(x)]);
  B1(i) = v.'*Bmat*v;
  B3(i) = v.'*Bmat*vb;
end
t1 = (2*sqrt(-1)*om*eye(numel(x))-A)\B1;
t2 = A\B3;
B2 = zeros(numel(x),1);
B4 = zeros(numel(x),1);
for i=1:numel(x)
  Bmat  = reshape(B(i,:,:),[numel(x),numel(x)]);
  B2(i) = vb.'*Bmat*t1;
  B4(i) = v.'*Bmat*t2;
end

y = real(w'*B2-2*w'*B4)/2/om;

end
