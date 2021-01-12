function y = piecewise(x, p, mode)
%PIECEWISE   'hspo'-compatible encoding of vector field.

x1 = x(1,:);
x2 = x(2,:);

r  = sqrt(x1.^2+x2.^2);

switch mode
  case 'left'
    y(1,:) = (1-r).*x1-x2;
    y(2,:) = x1+(1-r).*x2;
  case 'right'
    al = p(1,:);
    be = p(2,:);
    ga = p(3,:);
    
    y(1,:) = al.*(be-r).*x1-(ga+be-r).*x2;
    y(2,:) = al.*(be-r).*x2+(ga+be-r).*x1;
end

end
