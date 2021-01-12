function J = linode_DFDP(t, x, p) %#ok<INUSD,INUSL>
%LINODE_DFDP   'coll'-compatible encoding of Jacobian of 'linode' vector field w.r.t. parameters.
%
% Encoding is of an autonomous vector field.

x1 = x(1,:);

J = zeros(2,1,numel(x1));
J(2,1,:) = -x1;

end
