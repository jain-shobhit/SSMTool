function y = linode(x, p)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
om = p(1,:);

y(1,:) = x2;
y(2,:) = -x2 - x1 + cos(x3);
y(3,:) = om;

end
