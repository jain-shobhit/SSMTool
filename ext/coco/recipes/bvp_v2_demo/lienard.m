function y = lienard(x, p)
%LIENARD   'coll'-compatible encoding of lienard vector field.

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = p1.*x2-x2.^3-x1;

end
