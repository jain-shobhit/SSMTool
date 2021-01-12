function J = tor_DFDX(x, p)
%TOR_DFDX   'coll'-compatible encoding of Jacobian of vector field w.r.t. problem variables.

nu = p(1,:);
be = p(2,:);
ga = p(3,:);
r  = p(4,:);
a3 = p(5,:);
b3 = p(6,:);

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

J = zeros(3,3,numel(x1));

J(1,1,:) = ( -(be+nu) - 3*a3.*x1.^2 - 3*b3.*(x2-x1).^2 )./r;
J(1,2,:) = ( be + 3*b3.*(x2-x1).^2 )./r;
J(2,1,:) = be + 3*b3.*(x2-x1).^2;
J(2,2,:) = -(be+ga) - 3*b3.*(x2-x1).^2;
J(2,3,:) = -1;
J(3,2,:) = 1;

end
