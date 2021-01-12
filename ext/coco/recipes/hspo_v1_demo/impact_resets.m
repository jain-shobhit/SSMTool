function y = impact_resets(x, p, reset)
%IMPACT_RESETS   'hspo'-compatible encoding of reset function.

p6 = p(6,:);

y = x;
switch reset
  case 'bounce'
    y(2,:) = -p6.*y(2,:);
  case 'phase'
    y(3,:) = y(3,:)-2*pi;
end

end
