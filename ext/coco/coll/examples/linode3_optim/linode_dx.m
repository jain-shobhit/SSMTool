function J = linode_dx(x, p)

x3 = x(3,:);
k  = p(1,:);
th = p(2,:);

J = zeros(3,3,numel(k));
J(1,2,:) = 1;
J(2,1,:) = -k;
J(2,2,:) = -1;
J(2,3,:) = -sin(x3+th);

end
