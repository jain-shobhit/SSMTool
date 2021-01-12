function J = duffing_dp(x, p) %#ok<INUSD>

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

J = zeros(3,4,numel(x1));
J(2,1,:) = -2*x2;
J(2,2,:) = -x1.^3;
J(2,3,:) = cos(x3);
J(3,4,:) = 1;

end
