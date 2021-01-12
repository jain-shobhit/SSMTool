function dJ = linode_dxdp(t, x, p) %#ok<INUSD>

dJ = zeros(2,2,2,numel(t));
dJ(2,1,1,:) = -1;

end
