function J = linode_DFDP(x, p)
%LINODE_DFDP   'coll'-compatible encoding of Jacobian of 'linode' vector field w.r.t. parameters.

x1 = x(1,:);

J = zeros(3,1,numel(x1));
J(2,1,:) = -x1;

end
