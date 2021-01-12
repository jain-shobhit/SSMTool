function [data y] = lemniscate(prob, data, u)
%LEMNISCATE   COCO-compatible encoding of lemniscate equation.
  y = (u(1)^2+u(2)^2)^2+u(2)^2-u(1)^2;
end
