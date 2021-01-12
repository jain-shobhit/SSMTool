function dJ = duffing_dxdx(x, p)

x1 = x(1,:);
x3 = x(3,:);
p2 = p(2,:);
p3 = p(3,:);

dJ = zeros(3,3,3,numel(x1));
dJ(2,1,1,:) = -6*p2.*x1;
dJ(2,3,3,:) = -p3.*cos(x3);

end
