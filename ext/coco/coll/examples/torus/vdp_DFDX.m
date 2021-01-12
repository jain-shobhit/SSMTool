function J = vdp_DFDX(t, x, p)
%VDP_DFDX   'coll'-compatible encoding of Jacobian of langford vector field w.r.t. problem variables

x1 = x(1,:);
x2 = x(2,:);
c  = p(2,:);

J = zeros(2,2,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -1-2*c.*x1.*x2;
J(2,2,:) = c.*(1-x1.^2);

end
