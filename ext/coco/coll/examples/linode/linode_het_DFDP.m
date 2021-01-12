function J = linode_het_DFDP(t, x, p) %#ok<INUSD,INUSL>
%LINODE_HET_DFDP   'coll'-compatible encoding of Jacobian with respect to parameters
%
% Encoding is of an autonomous vector field.

x1 = x(1,:);

J = zeros(2,1,numel(x1));
J(2,1,:) = -x1;

end
