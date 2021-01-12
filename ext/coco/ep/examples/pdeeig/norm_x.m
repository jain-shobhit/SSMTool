function [data, y] = norm_x(opts, data, u) %#ok<INUSL>
%NORM_X   COCO-compatible function encoding of norm of input argument.
  y = norm(u);
end