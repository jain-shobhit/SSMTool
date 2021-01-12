function dJ = linode_dpdp(x, p)

x3 = x(3,:);
th = p(2,:);

dJ = zeros(3,2,2,numel(th));
dJ(2,2,2,:) = -cos(x3+th);

end
