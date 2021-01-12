function y = duff_events(x, p, event)
%DUFF_EVENTS   'hspo'-compatible encoding of duffing event function.

x3 = x(3,:);
om = p(5,:);

y = pi./om-x3; % Excitation switches every half period

end
