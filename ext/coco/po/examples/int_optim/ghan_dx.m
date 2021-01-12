function J = ghan_dx(x)

x1 = x(1,:);
x2 = x(2,:);

J  = zeros(1,3,numel(x1));
J(1,1,:) = 1./(1+x2.^2);
J(1,2,:) = -2*x1.*x2./(1+x2.^2).^2;

end
