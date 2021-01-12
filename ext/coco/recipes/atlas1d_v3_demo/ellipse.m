function [data y] = ellipse(prob, data, u)
%CIRCLE   COCO-compatible encoding of ellipse zero function.
  y = (2*(u(1)-1))^2+u(2)^2-1;
end
