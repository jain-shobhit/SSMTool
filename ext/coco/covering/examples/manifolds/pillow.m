function [data, y] = pillow(prob, data, u) %#ok<INUSL>
%PILLOW   COCO-compatible encoding of pillow zero function.
  y = u(1)^4 + u(2)^4 + u(3)^4 - u(1)^2 - u(2)^2 - u(3)^2;
end
