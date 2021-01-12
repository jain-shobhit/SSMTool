function y = impact_events(x, p, event)
%IMPACT_EVENTS   'hspo'-compatible encoding of event function.

x1 = x(1,:);
x3 = x(3,:);
p5 = p(5,:);

switch event
  case 'impact'
    y = p5-x1;
  case 'phase'
    y = pi-x3;
end

end
