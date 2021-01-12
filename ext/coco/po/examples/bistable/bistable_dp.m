function J = bistable_dp(t, x, p)
%BISTABLE_DP   'coll'-compatible encoding of Jacobian with respect to parameters.
%
% Encoding is of a non-autonomous vector field.

x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);

J = zeros(2,3,numel(p1));
J(2,1,:) = 2*pi./p1.^2.*p2.*t.*sin(2*pi./p1.*t);
J(2,2,:) = cos(2*pi./p1.*t);
J(2,3,:) = -x2;

end
