function J = cstr_dx(x, p)
%CSTR_DX Jacobian of vector field for CSTR model w.r.t. x
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
J = zeros(2,2,numel(z));

J(1,1,:) = -1 - l.*t.*z;
J(1,2,:) = l.*t.*(1-u).*z;
J(2,1,:) = -8*l.*t.*z;
J(2,2,:) = -1+8*l.*t.*(1-u).*z-t;

end
