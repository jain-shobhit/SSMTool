function [data, y] = doedel_evs(prob, data, u) %#ok<INUSL>
%DOEDEL_EVS   COCO-compatible encoding of doedel eigenspace conditions.

v = u(1:2);
w = u(3:4);
l = u(5);

y = [w-l*v; v'*v-1];

end
