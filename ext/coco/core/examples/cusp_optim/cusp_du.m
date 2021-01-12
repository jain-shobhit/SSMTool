function [data, J] = cusp_du(prob, data, u) %#ok<INUSL>
J = [3*u(1)^2-u(3), 1, -u(1)];
end
