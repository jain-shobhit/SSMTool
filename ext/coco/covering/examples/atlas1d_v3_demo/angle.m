function [data, y] = angle(prob, data, u) %#ok<INUSL>
%ANGLE   COCO-compatible encoding of angle zero function.

y = (u(1)-1)*sin(u(3))-u(2)*cos(u(3));

end
