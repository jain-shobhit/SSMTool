function [data, y] = ellipsoid(prob, data, u) %#ok<INUSL>
%ELLIPSOID   COCO-compatible encoding of ellipsoid zero function.
  y = (2*(u(1) - 1))^2 + u(2)^2 + u(3)^2 - 1;
end
