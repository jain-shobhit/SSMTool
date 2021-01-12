function [data y] = dist(prob, data, u)
%DIST   COCO-compatible encoding of squared distance to origin.
  y = u(1)^2+u(2)^2;
end
