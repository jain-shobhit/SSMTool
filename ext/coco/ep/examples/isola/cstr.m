function f = cstr(x, p)
%CSTR Vector field for CSTR model
%
% Continuous stirred tank reactor model of chemical reactions analyzed in
% "On the Numerical Continuation of Isolas of Equilibria," by Avitabile,
% Desroches, and Rodriquez, International Journal of Bifurcation and Chaos,
% 22(11), art. no. 1250277, 2012.

u = x(1,:);
v = x(2,:);
t = p(1,:);
l = p(2,:);

z = exp(v);

f = [-u + l.*t.*(1-u).*z; -v + 8*l.*t.*(1-u).*z-t.*v];

end
