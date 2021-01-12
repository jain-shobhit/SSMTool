function J = lorenz_DFDXDX(x, p)
%LORENZ_DFDXDX   'varcoll'-compatible encoding of first partials of vector-field Jacobian with respect to problem variables.

s = p(1,:);

J = zeros(3,3,3,numel(s));
J(2,3,1,:) = -1;
J(3,2,1,:) = 1;
J(3,1,2,:) = 1;
J(2,1,3,:) = -1;

end
