function Jt = bistable_dt(t, x, p)
%BISTABLE   'coll'-compatible encoding of Jacobian with respect to time
%
% Encoding is of a non-autonomous vector field.

p1 = p(1,:);
p2 = p(2,:);

Jt = zeros(2,numel(t));
Jt(2,:) = -2*pi*p2./p1.*sin(2*pi./p1.*t);

end
