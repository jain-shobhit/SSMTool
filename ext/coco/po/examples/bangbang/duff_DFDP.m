function J = duff_DFDP(x, p, mode)
%DUFF   'hspo'-compatible encoding of Jacobian w.r.t. problem parameters.

x1 = x(1,:);
x2 = x(2,:);

J = zeros(3,5,numel(x1));
J(2,1,:) = -x2;
J(2,2,:) = -x1;
J(2,3,:) = -x1.^3;
J(2,4,:) = 1;
switch mode
  case 'neg'
    J(2,4,:) = -1;
end

end
