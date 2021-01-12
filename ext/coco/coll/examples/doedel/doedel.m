function y = doedel(x, p)
%DOEDEL   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);

y(1,:) = 1-x1.^2;
y(2,:) = p1.*x1+p2.*x2;

end
