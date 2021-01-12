function [data, y] = linode_aut_bc(prob, data, u) %#ok<INUSL>
%LINODE_AUT_BC   'bvp'-compatible encoding of linode boundary conditions
%
% Encoding is of an autonomous vector field.

x0 = u(1:3);
x1 = u(4:6);

y = [x1(1:2)-x0(1:2); x1(3)-x0(3)-2*pi; x0(1)]; % periodicity in R2 x S1 on Poincare section

end
