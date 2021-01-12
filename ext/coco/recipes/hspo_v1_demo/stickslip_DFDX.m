function J = stickslip_DFDX(x, p, mode)
%STICKSLIP_DFDX   'hspo'-compatible encoding of Jacobian of vector field.

switch mode
  case 'stick'
    x1  = x(1,:);
    x3  = x(3,:);
    V   = p(1,:);
    c   = p(2,:);

    dF1 = 2*V.^2.*sin(x3).^2./(1-x1).^3-1;
    dF2 = -c;
    dF3 = 2*V.^2.*cos(x3).*sin(x3)./(1-x1).^2;
    
    J(:,:,:) = zeros(3,3,numel(x1));
    J(1,2,:) = 1;
    J(2,1,:) = dF1;
    J(2,2,:) = dF2;
    J(2,3,:) = dF3;
  case 'slip'
    x1  = x(1,:);
    x2  = x(2,:);
    x4  = x(4,:);
    V   = p(1,:);
    c   = p(2,:);
    
    dF2 = 2*V.^2.*sin(x4).^2./(1-x2).^3-1;
    dF3 = -c;
    dF4 = 2*V.^2.*cos(x4).*sin(x4)./(1-x2).^2;
    
    J(:,:,:) = zeros(4,4,numel(x1));
    J(2,3,:) = 1;
    J(1,2,:) = -dF2/5;
    J(1,3,:) = -dF3/5;
    J(1,4,:) = -dF4/5;
    J(3,2,:) = 6*dF2/5;
    J(3,3,:) = 6*dF3/5;
    J(3,4,:) = 6*dF4/5;
end

end
