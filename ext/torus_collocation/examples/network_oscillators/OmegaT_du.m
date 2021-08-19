function [data, J] = OmegaT_du(prob, data, u) %#ok<INUSL>
J = [u(2) u(1)];
end