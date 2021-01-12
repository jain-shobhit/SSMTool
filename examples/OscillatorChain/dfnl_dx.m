function Df = dfnl_dx(x)
n = length(x);
x = [0; x; 0];
g = 3*(x(2:end) - x(1:end-1)).^2;

D = [-g(2:end-1); 0]; 
E = g(1:end-1)+g(2:end);
F = [0; -g(2:end-1)]; 

Df = spdiags([D E F], -1:1, n,n);
