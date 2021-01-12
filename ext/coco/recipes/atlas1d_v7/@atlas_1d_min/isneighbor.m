function flag = isneighbor(atlas, chart1, chart2)
%ISNEIGHBOR   Check if two charts are neighbors.
%
% Neighboring charts have base points inside a cone centered on the other
% base point, aligned along the corresponding tangent vector, and with
% desired opening angle. In addition, the projection of each base point on
% to the other tangent vector lies within R from the other base point.
%
% Identical to atlas1d_v3.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: isneighbor.m 2839 2015-03-05 17:09:01Z fschild $

al  = atlas.cont.almax;
ta  = tan(al);
R   = atlas.cont.h;
x1  = chart1.x;
x2  = chart2.x;
dx  = x2-x1;
x1s = chart1.TS*(chart1.TS'*dx);
x2s = chart2.TS*(chart2.TS'*dx);
dst = [norm(x1s), norm(x2s), norm(dx-x1s), norm(dx-x2s), ...
       subspace(chart1.TS, chart2.TS)];
dstmx = [R, R, ta*norm(x1s), ta*norm(x2s), al];
flag  = all(dst<dstmx);

end
