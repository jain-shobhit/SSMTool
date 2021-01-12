function J = cusp_DFDP(x, p)
%CUSP_DFDP   'alg'-compatible encoding of Jacobian of cusp normal form w.r.t. problem parameters

x = x(1,:);

J = zeros(1,2,numel(x));
J(1,1,:) = 1;
J(1,2,:) = -x;

end
