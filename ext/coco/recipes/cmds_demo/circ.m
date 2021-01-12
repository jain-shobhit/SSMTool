function [data y] = circ(prob, data, u)
%CIRC   COCO-compatible encoding of circle function.
  y = u(1)^2+(u(2)-1)^2-1;
end
