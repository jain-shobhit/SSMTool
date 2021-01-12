function J = linode_het_DFDT(t, x, p)
%LINODE_HET_DFDT   'coll'-compatible encoding of Jacobian with respect to time
%
% Encoding is of a non-autonomous vector field.

x1 = x(1,:);

J = zeros(2,numel(x1));
J(2,:) = -sin(t);

end
