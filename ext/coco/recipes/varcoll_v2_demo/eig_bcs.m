function [data y] = eig_bcs(prob, data, u)
%EIG_BCS   COCO-compatible encoding of lorenz eigenspace conditions.

x10 = u(1:3);
x20 = u(4:6);
x30 = u(7:9);
s   = u(10);
r   = u(11);
vec = u(12:14);
eps = u(15:16);

evec = [(1-s+sqrt((1-s)^2+4*r*s))/2/r; 1; 0];
y = [x10-eps(1).*evec; x20-(x30+eps(2)*vec)];

end
