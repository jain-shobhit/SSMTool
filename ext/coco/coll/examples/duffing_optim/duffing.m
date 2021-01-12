function y = duffing(x, p)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
zeta  = p(1,:);
alpha = p(2,:);
amp   = p(3,:);
omega = p(4,:);

y(1,:) = x2;
y(2,:) = amp.*cos(x3)-2*zeta.*x2-x1-alpha.*x1.^3;
y(3,:) = omega;

end
