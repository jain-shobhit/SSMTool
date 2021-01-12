function y = lorenz(x, p)
%LORENZ   'coll'-compatible encoding of lorenz vector field.

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
s  = p(1,:);
r  = p(2,:);
b  = p(3,:);

y(1,:) = -s.*x1+s.*x2;
y(2,:) = -x1.*x3+r.*x1-x2;
y(3,:) = x1.*x2-b.*x3;

end
