function J = bykov_dp(x, p)
%BYKOV_DP Jacobian of vector field for BYK model w.r.t. p
%
% Bykov-Yablonskii-Kim model of oxidation of carbon monoxide on platinum
% analyzed in Tutorial IV: Two-parameter bifurcation analysis of equilibria
% and limit cycles with MATCONT, by Yu.A. Kuznetsov, September 20, 2011.
% (see http://www.staff.science.uu.nl/~kouzn101/NBA/LAB4.pdf)

p4 = p(4,:);
p7 = p(7,:);

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

z = 1 - x1 - x2 - x3;
J = zeros(3,7,numel(z));

J(1,1,:) = 2*z.^2;
J(1,3,:) = -x1.*x2;
J(1,5,:) = -2*x1.^2;
J(2,2,:) = z;
J(2,3,:) = -x1.*x2;
J(2,6,:) = -x2;
J(3,4,:) = z - p7.*x3;
J(3,7,:) = -p4.*x3;

end
