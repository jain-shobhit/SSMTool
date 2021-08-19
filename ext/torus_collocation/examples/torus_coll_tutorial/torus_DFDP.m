function J = torus_DFDP(t, x, p)
%TORUS_DFDP   'coll'-compatible encoding of Jacobian of vector field w.r.t. problem parameters.

x1 = x(1,:);
x2 = x(2,:);
om = p(1,:);

r = sqrt(x1.^2+x2.^2);
J = zeros(2,2,numel(x1));
J(1,1,:) = -r.*t.*x1.*sin(om.*t);
J(1,2,:) = -x2;
J(2,1,:) = -r.*t.*x2.*sin(om.*t);
J(2,2,:) = x1;

end
