function J = linode_het_DFDX(t, x, p) %#ok<INUSL>
%LINODE_HET_DFDX   'coll'-compatible encoding of Jacobian with respect to state
%
% Encoding is of a non-autonomous vector field.

x1 = x(1,:);
p1 = p(1,:);

J = zeros(2,2,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -p1;
J(2,2,:) = -1;

end
