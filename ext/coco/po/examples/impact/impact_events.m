function y = impact_events(x, p, event)
%IMPACT_EVENTS   'hspo'-compatible encoding of event function.

switch event
  case 'impact'
    y = p(5,:)-x(1,:);
  case 'phase'
    y = pi-x(3,:);
end

end
