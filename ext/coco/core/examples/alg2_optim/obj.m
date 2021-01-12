function [data, y] = obj(prob, data, u) %#ok<INUSL>
y = u(1)^2 + u(2);
end
