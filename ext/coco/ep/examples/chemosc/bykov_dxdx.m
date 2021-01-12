function J = bykov_dxdx(x, p) %#ok<INUSL>
%BYKOV_DXDX 2nd Jacobian of vector field for BYK model w.r.t. x
%
% Bykov-Yablonskii-Kim model of oxidation of carbon monoxide on platinum
% analyzed in Tutorial IV: Two-parameter bifurcation analysis of equilibria
% and limit cycles with MATCONT, by Yu.A. Kuznetsov, September 20, 2011.
% (see http://www.staff.science.uu.nl/~kouzn101/NBA/LAB4.pdf)

p1 = p(1,:);
p3 = p(3,:);
p5 = p(5,:);

J = zeros(3,3,3,numel(p1));

J(1,1,1,:) = 4*p1 - 4*p5;
J(1,1,2,:) = 4*p1 - p3;
J(1,1,3,:) = 4*p1;
J(1,2,1,:) = 4*p1 - p3;
J(1,2,2,:) = 4*p1;
J(1,2,3,:) = 4*p1;
J(1,3,1,:) = 4*p1;
J(1,3,2,:) = 4*p1;
J(1,3,3,:) = 4*p1;
J(2,1,2,:) = -p3;
J(2,2,1,:) = -p3;

end
