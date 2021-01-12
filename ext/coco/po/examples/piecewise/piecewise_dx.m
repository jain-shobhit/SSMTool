function J = piecewise_dx(x, p, mode)
%PIECEWISE_DFDX   'hspo'-compatible encoding of Jacobian of vector field w.r.t. state variables.

x1 = x(1,:);
x2 = x(2,:);

r  = sqrt(x1.^2+x2.^2);
rx = x1./r;
ry = x2./r;

J  = zeros(2,2,numel(r));
switch mode
  case 'left'
    J(1,1,:) = 1-r-x1.*rx;
    J(1,2,:) = -1-x1.*ry;
    J(2,1,:) = 1-x2.*rx;
    J(2,2,:) = 1-r- x2.*ry;
  case 'right'
    al   = p(1,:);
    be   = p(2,:);
    ga   = p(3,:);

    al_x = al.*x1-x2;
    al_y = al.*x2+x1;
    
    J(1,1,:) = al.*be-al.*r-rx.*al_x;
    J(1,2,:) = -be-ga+r-ry.*al_x;
    J(2,1,:) = be+ga-r-rx.*al_y;
    J(2,2,:) = al.*be-al.*r-ry.*al_y;
end

end
