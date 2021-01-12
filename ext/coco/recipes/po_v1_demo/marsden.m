function y = marsden(x, p)
%MARSDEN   'coll'-compatible encoding of marsden vector field.

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
p1 = p(1,:);
p2 = p(2,:);

y(1,:) = p1.*x1+x2+p2.*x1.^2;
y(2,:) = -x1+p1.*x2+x2.*x3;
y(3,:) = (p1.^2-1).*x2-x1-x3+x1.^2;

end
