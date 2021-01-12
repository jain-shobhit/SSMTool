function f = bykov(x, p)
%BYKOV Vector field for BYK model
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
f = [2*p1.*z.^2 - 2*p5.*x1.^2 - p3.*x1.*x2;
     p2.*z - p6.*x2 - p3.*x1.*x2;
     p4.*z - p7.*p4.*x3];

end
