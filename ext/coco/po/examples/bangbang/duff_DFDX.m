function J = duff_DFDX(x, p, mode) %#ok<INUSD>
%DUFF   'hspo'-compatible encoding of the Jacobian w.r.t. state variables.

x1 = x(1,:);
la = p(1,:);
al = p(2,:);
ep = p(3,:);

J = zeros(3,3,numel(x1));
J(1,2,:) = 1;
J(2,1,:) = -al-3*ep.*x1.^2;
J(2,2,:) = -la;

end
