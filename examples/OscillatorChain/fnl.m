function f = fnl(x)

x = [0; x; 0];
g = (x(2:end) - x(1:end-1)).^3; 

f = g(1:end-1)-g(2:end);