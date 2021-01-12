function y = torus(t, x, p)
%TORUS   'coll'-compatible encoding of vector field.

x1 = x(1,:);
x2 = x(2,:);
om = p(1,:);
fr = p(2,:);

r      = sqrt(x1.^2+x2.^2);
y(1,:) = -fr.*x2+x1.*(1+r.*(cos(om.*t)-1));
y(2,:) = fr.*x1+x2.*(1+r.*(cos(om.*t)-1));

end
