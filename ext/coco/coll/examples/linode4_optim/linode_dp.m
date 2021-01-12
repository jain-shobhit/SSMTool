function J = linode_dp(x, p) %#ok<INUSL>

om = p(1,:);

J = zeros(3,1,numel(om));
J(3,1,:) = 1;

end
