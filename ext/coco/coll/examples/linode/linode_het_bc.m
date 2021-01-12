function [data, y] = linode_het_bc(prob, data, u) %#ok<INUSL>
%LINODE_HET_BC   'bvp'-compatible encoding of linode boundary conditions
%
% Encoding is of a non-autonomous vector field.

x0 = u(1:2);
x1 = u(3:4);
T0 = u(5);
T  = u(6);

y = [x1(1:2)-x0(1:2); T0; T-2*pi]; % periodicity in R2 and fixed period

end
