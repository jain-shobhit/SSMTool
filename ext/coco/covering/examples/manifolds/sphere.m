function [data y] = sphere(prob, data, u) %#ok<INUSL>
%SPHERE   COCO-compatible encoding of sphere zero function.
  y = (u(1) - 1)^2 + u(2)^2 + u(3)^2 - 1;
end
