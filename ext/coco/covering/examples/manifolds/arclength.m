function [data, y] = arclength(prob, data, u) %#ok<INUSL>
%ARCLENGTH   COCO-compatible encoding of signed arclength zero function.
y = (u(1)-1)*sin(u(3))-u(2)*cos(u(3));
end
