function J = lorenz_DFDXDP(x, p)
%LORENZ_DFDXDP   'varcoll'-compatible encoding of first partials of vector-field Jacobian with respect to problem parameters.

s = p(1,:);

J = zeros(3,3,3,numel(s));
J(1,1,1,:) = -1;
J(1,2,1,:) = 1;
J(2,1,2,:) = 1;
J(3,3,3,:) = -1;

end
