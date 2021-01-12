function [data, f] = cusp(prob, data, u) %#ok<INUSL>
f = u(2)-u(1)*(u(3)-u(1)^2);
end
