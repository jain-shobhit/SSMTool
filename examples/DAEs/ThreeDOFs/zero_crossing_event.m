function [position,isterminal,direction] = zero_crossing_event(t,y)
  position = y(1); % The value that we want to be zero
  isterminal = 0;  % Halt integration 
  direction = 1;   % The zero can be approached from either direction
end