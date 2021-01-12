function J = tor_dp(x, p)
%TOR_DP   'coll'-compatible encoding of Jacobian of tor vector field.

x1 = x(1,:);
x2 = x(2,:);
nu = p(1,:);
be = p(2,:);
r  = p(4,:);
a3 = p(5,:);
b3 = p(6,:);

J = zeros(3,6,numel(x1));

J(1,1,:) = -x1./r;
J(1,2,:) = (-x1+x2)./r;
J(1,4,:) = -( -(be+nu).*x1 + be.*x2 - a3.*x1.^3 + b3.*(x2-x1).^3 )./r.^2;
J(1,5,:) = -x1.^3./r;
J(1,6,:) = (x2-x1).^3./r;
J(2,2,:) = x1-x2;
J(2,3,:) = -x2;
J(2,6,:) = -(x2-x1).^3;

end
