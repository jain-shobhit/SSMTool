function [data y] = angle(prob, data, u)
%ANGLE   COCO-compatible encoding of angle monitor function.
  y = atan(u(3)./(u(1)+0.2));
end
