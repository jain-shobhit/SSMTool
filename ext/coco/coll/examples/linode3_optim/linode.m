function y = linode(x, p)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
k  = p(1,:);
th = p(2,:);

y(1,:) = x2;
y(2,:) = -x2-k.*x1+cos(x3+th);
y(3,:) = 1;

end
