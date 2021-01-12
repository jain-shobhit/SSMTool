function dJ = linode_dxdp(x, p)

x3 = x(3,:);
th = p(2,:);

dJ = zeros(3,3,2,numel(th));
dJ(2,1,1,:) = -1;
dJ(2,3,2,:) = -cos(x3+th);

end
