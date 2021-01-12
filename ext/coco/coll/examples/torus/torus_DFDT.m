function J = torus_DFDT(t, x, p)
%TORUS_DFDT   'coll'-compatible encoding of vector field.

x1 = x(1,:);
x2 = x(2,:);
om = p(1,:);

r = sqrt(x1.^2+x2.^2);
J = zeros(2,numel(x1));
J(1,:) = -x1.*r.*om.*sin(om.*t);
J(2,:) = -x2.*r.*om.*sin(om.*t);

end
