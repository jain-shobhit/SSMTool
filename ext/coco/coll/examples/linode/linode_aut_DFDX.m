function J = linode_aut_DFDX(x, p)
%LINODE_AUT_DFDX   'coll'-compatible encoding of Jacobian with respect to state
%
% Encoding is of an autonomous vector field.

x1 = x(1,:);
x3 = x(3,:);
p1 = p(1,:);

J = zeros(3,3,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -p1;
J(2,2,:) = -1;
J(2,3,:) = -sin(x3);

end
