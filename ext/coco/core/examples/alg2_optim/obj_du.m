function [data, J] = obj_du(prob, data, u) %#ok<INUSL>
J = [2*u(1) 1];
end
