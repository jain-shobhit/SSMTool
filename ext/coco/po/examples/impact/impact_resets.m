function y = impact_resets(x, p, reset)
%IMPACT_RESETS   'hspo'-compatible encoding of reset function.

y = x;

switch reset
  case 'bounce'
    y(2,:) = -p(6,:).*y(2,:);
  case 'phase'
    y(3,:) = y(3,:)-2*pi;
end

end
