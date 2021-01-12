function [data y] = empty(prob, data, u)
%EMPTY   COCO-compatible encoding of everywhere positive zero function.
  y = u(1)^2+u(2)^2+1;
end
