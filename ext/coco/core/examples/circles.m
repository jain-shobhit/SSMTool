function [data y] = circles(prob, data, u)
y = (u(1)^2+u(2)^2-1)*((u(1)-1)^2+u(2)^2-1);
end
