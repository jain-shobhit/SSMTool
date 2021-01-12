function [data, f] = sphere(prob, data, u) %#ok<INUSL>
f = u(1)^2 + u(2)^2 + u(3)^2 + u(4)^2 - 1;
end
