function [data, y] = linode_bc(prob, data, u) %#ok<INUSL>

x0 = u(1:3);
x1 = u(4:6);

y = [x1(1:2)-x0(1:2); x1(3)-x0(3)-2*pi; x0(3)];

end
