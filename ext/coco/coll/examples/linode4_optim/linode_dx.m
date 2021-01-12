function J = linode_dx(x, p) %#ok<INUSD>

x3 = x(3,:);

J = zeros(3,3,numel(x3));
J(1,2,:) = 1;
J(2,1,:) = -1;
J(2,2,:) = -1;
J(2,3,:) = -sin(x3);

end
