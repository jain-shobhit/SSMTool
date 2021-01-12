function y = bistable(t, x, p)
%BISTABLE   'coll'-compatible encoding of harmonically excited hardening oscillator
%
% Encoding is of a non-autonomous vector field.

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);

y(1,:) = x2;
y(2,:) = -p3.*x2-x1-x1.^3+p2.*cos(2*pi./p1.*t);

end
