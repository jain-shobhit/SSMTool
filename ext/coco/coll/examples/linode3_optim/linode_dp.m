function J = linode_dp(x, p)

x1 = x(1,:);
x3 = x(3,:);
th = p(2,:);

J = zeros(3,2,numel(th));
J(2,1,:) = -x1;
J(2,2,:) = -sin(x3+th);

end
