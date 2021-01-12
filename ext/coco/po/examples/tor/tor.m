function y = tor(x, p)
%TOR   'coll'-compatible encoding of tor vector field.

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
nu = p(1,:);
be = p(2,:);
ga = p(3,:);
r  = p(4,:);
a3 = p(5,:);
b3 = p(6,:);

y(1,:) = ( -(be+nu).*x1 + be.*x2 - a3.*x1.^3 + b3.*(x2-x1).^3 )./r;
y(2,:) =  be.*x1 - (be+ga).*x2 - x3 - b3.*(x2-x1).^3;
y(3,:) = x2;

end
