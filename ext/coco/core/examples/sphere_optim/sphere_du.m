function [data, J] = sphere_du(prob, data, u) %#ok<INUSL>
J = [2*u(1), 2*u(2), 2*u(3), 2*u(4)];
end
