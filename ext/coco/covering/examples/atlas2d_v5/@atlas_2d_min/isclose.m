function flag = isclose(atlas, chart1, chart2)
%ISCLOSE   Check if two charts are close to each other.
%
% Symmetric test for closeness that generalizes isneighbor to consider
% charts that could have been obtained by two consecutive steps of
% continuation.
%
% Identical to atlas2d_v4.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: isclose.m 3087 2019-04-04 19:54:09Z hdankowicz $

al   = atlas.cont.almax;
R    = atlas.cont.R;
ta   = tan(al);
t2a  = tan(2*al);
x1   = chart1.xp;
x2   = chart2.xp;
dx   = x2-x1;
phi1 = chart1.TSp'*dx;
phi2 = chart2.TSp'*dx;
x1s  = chart1.TSp*(phi1);
x2s  = chart2.TSp*(phi2);
dst  = [norm(x1s), norm(x2s), norm(dx-x1s), norm(dx-x2s), ...
        subspace(chart1.TSp, chart2.TSp)];
n1mx  = ta*min(R,norm(x1s))+t2a*max(0,norm(x1s)-R);
n2mx  = ta*min(R,norm(x2s))+t2a*max(0,norm(x2s)-R);
dstmx = [2*R, 2*R, n1mx, n2mx, 2*al];
flag  = false;
if all(dst<dstmx)
  test1 = chart1.v.*(chart1.s'*phi1)-norm(phi1)^2/2;
  test2 = chart2.v.*(chart2.s'*phi2)+norm(phi2)^2/2;
  flag  = any(test1>0) && any(test2<0);
end

end
