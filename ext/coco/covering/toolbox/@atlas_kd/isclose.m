function flag = isclose(atlas, chart1, chart2)
%ISCLOSE   Check if two charts are close to each other.
%
% Symmetric test for closeness that generalizes isneighbor to consider
% charts that could have been obtained by two consecutive steps of
% continuation.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

al   = atlas.cont.almax;
R1   = chart1.R;
R2   = chart2.R;
ta   = tan(al);
t2a  = tan(2*al);
x1   = chart1.xp;
x2   = chart2.xp;
dx   = x2 - x1;
phi1 = chart1.TSp'*dx;
phi2 = chart2.TSp'*dx;
x1s  = chart1.TSp*(phi1);
x2s  = chart2.TSp*(phi2);
dst  = [norm(x1s), norm(x2s), norm(dx - x1s), norm(dx - x2s), ...
        subspace(chart1.TSp, chart2.TSp)];
n1mx  = ta*min(R1,norm(x1s)) + t2a*max(0,norm(x1s) - R1);
n2mx  = ta*min(R2,norm(x2s)) + t2a*max(0,norm(x2s) - R2);
dstmx = [R1+R2, R1+R2, n1mx, n2mx, 2*al];
flag  = all(dst<dstmx); %false;%
% if all(dst<dstmx)
%   test1 = cell2mat(chart1.P.v)*(chart1.s'*phi1)-norm(phi1)^2/2;
%   test2 = cell2mat(chart2.P.v)*(chart2.s'*phi2)+norm(phi2)^2/2;
%   flag  = any(test1>0) && any(test2<0);
% end

end
