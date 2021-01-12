function J = lang_DFDP(x, p)

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

J = zeros(3,3,numel(x1));
J(1,1,:) = -x2;
J(2,1,:) = x1;
J(3,2,:) = -x3.*(x1.^2+x2.^2);
J(3,3,:) = x3.*x1.^3;

end