function [data, y] = wdist(prob, data, u) %#ok<INUSL>
%WDIST Arclength and angle condition along isola
%
% Impose weighted arclength conditions on each pair of points along isola.

pt = reshape(u(1:end-1), data.shp);
ds = u(end)*ones(data.np,1);

dw = (pt - pt(data.shf)).*data.w;
y  = sqrt(sum(dw.^2, 1))' - ds;

end
