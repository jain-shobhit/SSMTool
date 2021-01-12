function [data, dJ] = sphere_dudu(prob, data, u) %#ok<INUSD,INUSL>

dJ = zeros(1,4,4);
dJ(1,1,1) = 2;
dJ(1,2,2) = 2;
dJ(1,3,3) = 2;
dJ(1,4,4) = 2;

end
