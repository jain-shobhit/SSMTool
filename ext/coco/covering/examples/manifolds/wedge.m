function [data, y] = wedge(prob, data, u) %#ok<INUSL>
%ANGLE   COCO-compatible encoding of angle monitor function.
  y = atan(u(3)./(u(1)+0.2));
end
