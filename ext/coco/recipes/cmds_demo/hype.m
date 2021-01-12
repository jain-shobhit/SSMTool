function [data y] = hype(prob, data, u)
%HYPE   COCO-compatible encoding of difference of squares.
  y = u(1)^2-u(2)^2;
end
