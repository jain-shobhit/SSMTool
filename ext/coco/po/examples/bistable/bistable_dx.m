function J = bistable_dx(t, x, p) %#ok<INUSL>
%NONLINODE_DFDX   'coll'-compatible encoding of Jacobian with respect to state.
%
% Encoding is of a non-autonomous vector field.

x1 = x(1,:);
p3 = p(3,:);

J = zeros(2,2,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -1-3*x1.^2;
J(2,2,:) = -p3;

end
