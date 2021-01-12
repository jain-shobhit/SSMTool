function J = duffing_dx(x, p)

x1 = x(1,:);
x3 = x(3,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);

J = zeros(3,3,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -1-3*p2.*x1.^2;
J(2,2,:) = -2*p1;
J(2,3,:) = -p3.*sin(x3);

end
