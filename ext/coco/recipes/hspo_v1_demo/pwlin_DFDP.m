function J = pwlin_DFDP(x, p, mode)
%PWLIN_DFDP   'hspo'-compatible encoding of Jacobian of vector field w.r.t. problem parameters.

x1 = x(1,:);
x2 = x(2,:);

r  = sqrt(x1.^2+x2.^2);

J  = zeros(2,3,numel(r));
switch mode
  case 'right'
    al = p(1,:);
    be = p(2,:);
    
    J(1,1,:) = (be-r).*x1;
    J(1,2,:) = al.*x1-x2;
    J(1,3,:) = -x2;
    J(2,1,:) = (be-r).*x2;
    J(2,2,:) = al.*x2+x1;
    J(2,3,:) = x1;
end

end
