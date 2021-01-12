function J = linode_aut_DFDP(x, p) %#ok<INUSD>
%LINODE_AUT_DFDP   'coll'-compatible encoding of Jacobian with respect to parameters
%
% Encoding is of an autonomous vector field.

x1 = x(1,:);

J = zeros(3,1,numel(x1));
J(2,1,:) = -x1;

end
