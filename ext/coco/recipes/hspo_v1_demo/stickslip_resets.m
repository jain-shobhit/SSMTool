function y = stickslip_resets(x, p, reset)
%STICKSLIP_RESETS   'hspo'-compatible encoding of reset function.

switch reset % corrected typo in Recipes for Continuation, 1st edition, page 240
  case 'bounce'
    e = p(5);
    
    y(1,:) = (1+e)/6*x(2,:);
    y(2,:) = x(1,:);
    y(3,:) = -e*x(2,:);
    y(4,:) = x(3,:);
  case 'phase'
    y(1,:) = x(1,:);
    y(2,:) = x(2,:);
    y(3,:) = x(3,:)-pi;
  case 'turn'
    y = x;
  case 'stick'
    y = x(2:4,:);
end

end
