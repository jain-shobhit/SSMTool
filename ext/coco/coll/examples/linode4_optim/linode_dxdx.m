function dJ = linode_dxdx(x, p) %#ok<INUSD>

x3 = x(3,:);

dJ = zeros(3,3,3,numel(x3));
dJ(2,3,3,:) = -cos(x3);

end
