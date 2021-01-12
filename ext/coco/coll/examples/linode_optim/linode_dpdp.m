function dJ = linode_dpdp(t, x, p) %#ok<INUSL>

th = p(2,:);

dJ = zeros(2,2,2,numel(t));
dJ(2,2,2,:) = -cos(t+th);

end
