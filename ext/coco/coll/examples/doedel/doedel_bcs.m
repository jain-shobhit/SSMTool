function [data, y] = doedel_bcs(prob, data, u) %#ok<INUSL>
%DOEDEL_BCS   COCO-compatible encoding of doedel boundary conditions.

x10 = u(1:2);
x20 = u(3:4);
eqs = u(5:8);
vec = u(9:10);
eps = u(11:12);
th  = u(13);

y = [x10-(eqs(1:2)+eps(1)*[cos(th); sin(th)]);...
     x20-(eqs(3:4)+eps(2)*vec)];

end
