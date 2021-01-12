function dJ = duffing_dxdp(x, p) %#ok<INUSD>

x1 = x(1,:);
x3 = x(3,:);

dJ = zeros(3,3,4,numel(x1));
dJ(2,1,2,:) = -3*x1.^2;
dJ(2,2,1,:) = -2;
dJ(2,3,3,:) = -sin(x3);

end
