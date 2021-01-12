function y = linode(t, x, p)
%LINODE   'coll'-compatible encoding of linode vector field.
%
% Encoding is of a non-autonomous vector field.

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);

y(1,:) = x2;
y(2,:) = -x2-p1.*x1+p2.*cos(t);

end
