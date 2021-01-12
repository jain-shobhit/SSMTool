function J = linode_dp(t, x, p)

x1 = x(1,:);
th = p(2,:);

J = zeros(2,2,numel(t));
J(2,1,:) = -x1;
J(2,2,:) = -sin(t+th);

end
