function [data y] = doedel_evs(prob, data, u)
%DOEDEL_EVS   COCO-compatible encoding of doedel eigenspace conditions.

par = u(1:2);
eqs = u(3:4);
vec = u(5:6);
lam = u(7);

jac = doedel_DFDX(eqs, par);

y = [jac*vec-lam*vec; vec'*vec-1];

end
