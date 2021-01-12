function [prob, atlas, cseg, flush] = init_atlas(atlas, prob, cseg)
%INIT_ATLAS   Initialize atlas.
%
% Assign successfully located initial chart to point list and flush.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

chart          = cseg.curr_chart; % this should probably be cseg.ptlist{1}
R              = atlas.cont.R;

% Initialize chart radius to keep residuum below MaxRes
[prob, chart, res] = chart_hc_res(prob, chart, R/atlas.cont.Rmarg);
while (res > atlas.cont.MaxRes) && (R > atlas.cont.R_min)
  R   = max(atlas.cont.R_min, atlas.cont.R_fac_min*R);
  [prob, chart, res] = chart_hc_res(prob, chart, R/atlas.cont.Rmarg);
end

chart.R        = R;
chart.h        = R;
mark           = chart.P.mark;
chart          = create_polyhedron(chart, atlas.cont.Rmarg);
chart.P.mark   = mark;
chart.bv       = find(~chart.P.mark);
chart.loc      = [];
chart.intsct   = [];
chart.pt_type  = 'EP';
cseg.ptlist{1} = chart;
flush          = true;

end

function [prob, chart, res] = chart_hc_res(prob, chart, h)

% Ignore residuum of arc-length constraint
prob = coco_emit(prob, 'set_mode', 'check_res');

% Compute residuum on boundary box
xx  = chart.TS*chart.G*chart.s(:,chart.bv);
res = 0;
for i=1:size(xx,2)
  [prob, chart, f] = prob.efunc.F(prob, chart, chart.x + h*xx(:,i));
  res = max(res,norm(f));
end

end
