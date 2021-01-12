function J = doedel_DFDX(x, p)
%DOEDEL_DFDX   'coll'-compatible encoding of Jacobian with respect to state

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);

J = zeros(2,2,numel(x1));
J(1,1,:) = -2*x1;
J(2,1,:) = p1;
J(2,2,:) = p2;

end
