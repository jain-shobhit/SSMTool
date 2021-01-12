function J = doedel_DFDP(x, p) %#ok<INUSD>
%DOEDEL_DFDP   'coll'-compatible encoding of Jacobian with respect to parameters

x1 = x(1,:);
x2 = x(2,:);

J = zeros(2,2,numel(x1));
J(2,1,:) = x1;
J(2,2,:) = x2;

end
