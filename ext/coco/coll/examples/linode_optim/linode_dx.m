function J = linode_dx(t, x, p) %#ok<INUSL>

k = p(1,:);

J = zeros(2,2,numel(t));
J(1,2,:) = 1;
J(2,1,:) = -k;
J(2,2,:) = -1;

end
