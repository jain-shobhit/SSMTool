function J = vdp_DFDP(t, x, p)
%VDP_DFDP   'coll'-compatible encoding of Jacobian of langford vector field w.r.t. problem parameters.

x1 = x(1,:);
x2 = x(2,:);
om = p(1,:);
a  = p(3,:);

J = zeros(2,3,numel(x1));
J(2,1,:) = -a.*t.*sin(om.*t);
J(2,2,:) = x2.*(1-x1.^2);
J(2,3,:) = cos(om.*t);

end
