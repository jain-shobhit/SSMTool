function J = cusp_DFDX(x, p)
%CUSP_DFDX   'alg'-compatible encoding of Jacobian of cusp normal form w.r.t. problem variables

x  = x(1,:);
la = p(2,:);

J = 3*x.^2-la;

end
