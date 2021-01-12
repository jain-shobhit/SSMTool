function [data, dJ] = cusp_dudu(prob, data, u) %#ok<INUSL>

dJ = zeros(1,3,3);
dJ(1,1,1) = 6*u(1);
dJ(1,1,3) = -1;
dJ(1,3,1) = -1;

end
