function flag = isneighbor(atlas, chart1, chart2)
%ISNEIGHBOR   Check if two charts are neighbors.
%
% Neighboring charts have base points inside a cone centered on the other
% base point, aligned along the corresponding tangent vector, and with
% desired opening angle. In addition, the projection of each base point on
% to the other tangent vector lies within R from the other base point.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: isneighbor.m 3087 2019-04-04 19:54:09Z hdankowicz $

al  = atlas.cont.almax;
ta  = tan(al);
R   = atlas.cont.R;
x1  = chart1.xp;
x2  = chart2.xp;
dx  = x2-x1;
x1s = chart1.TSp*(chart1.TSp'*dx);
x2s = chart2.TSp*(chart2.TSp'*dx);
dst = [norm(x1s), norm(x2s), norm(dx-x1s), norm(dx-x2s), ...
       subspace(chart1.TSp, chart2.TSp)];
dstmx = [R, R, ta*norm(x1s), ta*norm(x2s), al];
flag  = all(dst<dstmx);

end
