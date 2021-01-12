function [data, y] = angles(prob, data, u) %#ok<INUSL>
%ANGLE   COCO-compatible encoding of angles zero function.

y = [u(1)*sin(u(4))-u(2)*cos(u(4)); ...
  (sqrt(u(1)^2+u(2)^2)-2)*sin(u(5))-u(3)*cos(u(5))];

end
