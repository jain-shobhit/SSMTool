function [data, y] = linode_bc(prob, data, u) %#ok<INUSL>

x0 = u(1:2);
x1 = u(3:4);
T0 = u(5);
T  = u(6);

y = [x1(1:2)-x0(1:2); T0; T-2*pi];

end
