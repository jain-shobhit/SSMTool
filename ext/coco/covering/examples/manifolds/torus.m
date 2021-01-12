function [data, y] = torus(prob, data, u) %#ok<INUSL>
%TORUS   COCO-compatible encoding of torus zero function.
  y = (2 - sqrt(u(1)^2 + u(2)^2))^2 + u(3)^2 - 1;
end
