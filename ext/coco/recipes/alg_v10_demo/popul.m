function y = popul(x, p)
%POPUL   'alg'-compatible encoding of population dynamics model
  y = [x(1)-x(1)*x(2)/(1+p(1)*x(1)); x(1)*x(2)/(1+p(1)*x(1))-x(2)-p(2)*x(2)^2];
end
