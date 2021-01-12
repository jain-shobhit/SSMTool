function [data y] = euclid(prob, data, u)
%EUCLID   COCO-compatible encoding of distance to origin.
  y = sqrt(u(1)^2+u(2)^2);
end
