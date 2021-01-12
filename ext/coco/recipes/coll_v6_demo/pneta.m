function y = pneta(x, p)
%PNETA   'coll'-compatible encoding of pneta vector field.

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);

y(1,:) = x2;
y(2,:) = 0.5*p1.*x2-p1.*x2.^3-x1;

end
