function dJ = linode_dxdx(x, p)

x3 = x(3,:);
th = p(2,:);

dJ = zeros(3,3,3,numel(th));
dJ(2,3,3,:) = -cos(x3+th);

end
