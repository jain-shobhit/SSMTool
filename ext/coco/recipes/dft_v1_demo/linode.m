function y = linode(xx, pp)
%LINODE   'dft'-compatible encoding of vector field.

p1 = pp(1,:);

x1 = xx(1,:);
x2 = xx(2,:);
x3 = xx(3,:);
x4 = xx(4,:);

y(1,:) = x2;
y(2,:) = -x2 - p1.*x1 + x4;
y(3,:) =  x3 + x4 - x3.*(x3.^2+x4.^2);
y(4,:) =  x4 - x3 - x4.*(x3.^2+x4.^2);

end
