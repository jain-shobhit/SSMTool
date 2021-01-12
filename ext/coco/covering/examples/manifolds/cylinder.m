function [data y] = cylinder(prob, data, u) %#ok<INUSL>
%CYLINDER   COCO-compatible encoding of cylinder zero function.
  y = (u(1)-1)^2 + u(2)^2 - 1;
end
