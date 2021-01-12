function y = hopf(x, p)
%HOPF   'coll'-compatible encoding of hopf vector field.

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);

r2 = x1.^2+x2.^2;

y(1,:) = x1.*(p1+p2.*r2-r2.^2)-x2;
y(2,:) = x2.*(p1+p2.*r2-r2.^2)+x1;

end
