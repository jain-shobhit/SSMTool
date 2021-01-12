function J = bykov_dx(x, p)
%BYKOV_DX Jacobian of vector field for BYK model w.r.t. x
%
% Bykov-Yablonskii-Kim model of oxidation of carbon monoxide on platinum
% analyzed in Tutorial IV: Two-parameter bifurcation analysis of equilibria
% and limit cycles with MATCONT, by Yu.A. Kuznetsov, September 20, 2011.
% (see http://www.staff.science.uu.nl/~kouzn101/NBA/LAB4.pdf)

p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
p4 = p(4,:);
p5 = p(5,:);
p6 = p(6,:);
p7 = p(7,:);

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);

z = 1 - x1 - x2 - x3;
J = zeros(3,3,numel(z));

J(1,1,:) = -4*p5.*x1 - p3.*x2 - 4*p1.*z;
J(1,2,:) = -p3.*x1 - 4*p1.*z;
J(1,3,:) = -4*p1.*z;
J(2,1,:) = -p2 - p3.*x2;
J(2,2,:) = -p2 - p6 - p3.*x1;
J(2,3,:) = -p2;
J(3,1,:) = -p4;
J(3,2,:) = -p4;
J(3,3,:) = -p4 - p4.*p7;

end
