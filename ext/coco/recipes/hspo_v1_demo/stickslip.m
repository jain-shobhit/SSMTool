function y = stickslip(x, p, mode)
%STICKSLIP   'hspo'-compatible encoding of vector field.

switch mode
  case 'stick'
    x1 = x(1,:);
    x2 = x(2,:);
    x3 = x(3,:);
    V  = p(1,:);
    c  = p(2,:);
    w  = p(4,:);

    F  = V.^2.*sin(x3).^2./(1-x1).^2-c.*x2-x1;
    
    y(1,:) = x2;
    y(2,:) = F;
    y(3,:) = w;
  case 'slip'
    x2 = x(2,:);
    x3 = x(3,:);
    x4 = x(4,:);
    V  = p(1,:);
    c  = p(2,:);
    n  = p(3,:);
    w  = p(4,:);
    
    F  = V.^2.*sin(x4).^2./(1-x2).^2-c.*x3-x2;
    
    y(1,:) = -n-F/5;
    y(2,:) = x3;
    y(3,:) = n+6*F/5;
    y(4,:) = w;
end

end
