function chart1 = chop_polyhedron(edge, chart1, chart2)
% CHOP_POLYHEDRON Modify the polyhedron
%
% CHART1 = CHOP_POLYHEDRON(EDGE, CHART1, CHART2)
%
% EDGE   - Edge number
% CHART1 - Chart to be modified
% CHART2 - Intersecting chart
%
% Modify the polyhedron associated with chart1 to reflect the intersection
% of its defining hypersphere with that of chart2.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

R1   = chart1.R;
R2   = chart2.R;
dir  = chart2.xp - chart1.xp;
nrm  = chart1.TSp'*dir;
dist = norm(nrm);
nrm  = nrm/dist;
onrm = .5*(dist*dist + R1*R1 - R2*R2)/dist;

% Cut polyhedron by half space
chart1.nrm        = [chart1.nrm      nrm   ];
chart1.onrm       = [chart1.onrm    onrm   ];
chart1.neigh      = [chart1.neigh chart2.id];
[chart1.P, rem]   = subtract_half_space(chart1.P, edge, nrm, onrm);
chart1.nrm(:,rem) = [];
chart1.onrm(rem)  = [];
chart1.neigh(rem) = [];

% Update index array of available vertices
idx = cellfun(@(v) norm(v), chart1.P.v) < chart1.R;
chart1.P.mark(idx) = 1;
chart1.bv = find(~chart1.P.mark);
chart1.s  = cell2mat(cellfun(@(v) v/norm(v), chart1.P.v, ...
  'UniformOutput', false));

end
