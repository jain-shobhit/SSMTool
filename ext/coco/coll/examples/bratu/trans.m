function [data, y] = trans(prob, data, u) %#ok<INUSL>
%TRANS Associate the problem parameter p with the integration constant C

y = u(1)-4*u(2).^2/(1+cosh(u(2)));

end
