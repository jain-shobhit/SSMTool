function J = linode_dt(t, x, p)  %#ok<INUSL>

th = p(2,:);

J = zeros(2,numel(t));
J(2,:) = -sin(t+th);

end
