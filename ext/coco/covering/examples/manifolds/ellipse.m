function [data, y] = ellipse(prob, data, u) %#ok<INUSL>
%ELLIPSE   COCO-compatible encoding of ellipse zero function.
  y = (10*(u(1) - 1))^2 + u(2)^2 - 1;
end
