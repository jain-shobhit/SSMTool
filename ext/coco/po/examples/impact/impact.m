function y = impact(x, p, mode)
%IMPACT   'hspo'-compatible encoding of vector field.

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
p4 = p(4,:);

y(1,:) = x2;
y(2,:) = p3.*cos(x3)-p2.*x2-p1.*x1;
y(3,:) = p4;

end
