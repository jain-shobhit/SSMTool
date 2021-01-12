function J = lorenz_DFDP(x, p)
%LORENZ_DFDP   'coll'-compatible encoding of Jacobian of lorenz vector field w.r.t. parameters.

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

J = zeros(3,3,numel(x1));
J(1,1,:) = -x1+x2;
J(2,2,:) = x1;
J(3,3,:) = -x3;

end
