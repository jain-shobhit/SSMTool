function y = stickslip_events(x, p, event) %#ok<INUSL>
%STICKSLIP_EVENTS   'hspo'-compatible encoding of event function.

switch event % corrected typo in Recipes for Continuation, 1st edition, page 240
  case 'collision'
    y = 0.5-x(1,:);
  case 'phase'
    y = pi/2-x(3,:);
  case 'minsep'
    y = x(2,:);
  case 'rest'
    y = x(1,:);
end

end
