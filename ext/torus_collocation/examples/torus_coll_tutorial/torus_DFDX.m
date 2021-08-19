function J = torus_DFDX(t, x, p)
%TORUS_DFDX   'coll'-compatible encoding of Jacobian of vector field w.r.t. problem variables

x1 = x(1,:);
x2 = x(2,:);
om = p(1,:);
fr = p(2,:);

r = sqrt(x1.^2+x2.^2);
J = zeros(2,2,numel(x1));
J(1,1,:) = (r+(r.^2+x1.^2).*(cos(om.*t)-1))./r;
J(1,2,:) = (-fr.*r+x1.*x2.*(cos(om.*t)-1))./r;
J(2,1,:) = (fr.*r+x1.*x2.*(cos(om.*t)-1))./r;
J(2,2,:) = (r+(r.^2+x2.^2).*(cos(om.*t)-1))./r;

end
