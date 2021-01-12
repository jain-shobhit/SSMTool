function [x0, y0] = compute_fixed_points_2D(x,y,xd,yd)
% This function returns the fixed points of a 2D vector field [xd, yd]  
% of in the mesh grid defined by x, y. The output [x0, y0] are the set of 
% points where (xd, yd) is numerically zero.

% zero level set of xd 
r1 = contourc(x,y,xd,[0,0]);
[x10,y10,~] = get_contour_xy(r1);

% zero level set of yd 
r2 = contourc(x,y,yd,[0,0]);
[x20,y20,~] = get_contour_xy(r2);

% intersection of two level curves
[x0, y0] = polyxpoly(x10,y10,x20,y20);
end