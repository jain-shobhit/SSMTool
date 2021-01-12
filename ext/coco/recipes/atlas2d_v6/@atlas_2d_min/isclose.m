function flag = isclose(atlas, chart1, chart2)
%ISCLOSE   Check if two charts are close to each other.
%
% Symmetric test for closeness that generalizes isneighbor to consider
% charts that could have been obtained by two consecutive steps of
% continuation.
%
% Identical to atlas2d_v4.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: isclose.m 2839 2015-03-05 17:09:01Z fschild $

al   = atlas.cont.almax;
R    = atlas.cont.h;
ta   = tan(al);
t2a  = tan(2*al);
x1   = chart1.x;
x2   = chart2.x;
dx   = x2-x1;
phi1 = chart1.TS'*dx;
phi2 = chart2.TS'*dx;
x1s  = chart1.TS*(phi1);
x2s  = chart2.TS*(phi2);
dst  = [norm(x1s), norm(x2s), norm(dx-x1s), norm(dx-x2s), ...
        subspace(chart1.TS, chart2.TS)];
n1mx  = ta*min(R,norm(x1s))+t2a*max(0,norm(x1s)-R);
n2mx  = ta*min(R,norm(x2s))+t2a*max(0,norm(x2s)-R);
dstmx = [2*R, 2*R, n1mx, n2mx, 2*al];
flag  = false;
if all(dst<dstmx);
  test1 = chart1.v.*(chart1.s'*phi1)-norm(phi1)^2/2;
  test2 = chart2.v.*(chart2.s'*phi2)+norm(phi2)^2/2;
  flag  = any(test1>0) && any(test2<0);
end

end
