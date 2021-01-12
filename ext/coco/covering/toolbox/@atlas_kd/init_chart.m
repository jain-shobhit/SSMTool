function [prob, atlas, cseg] = init_chart(atlas, prob, cseg)
%INIT_CHART   Initialize chart.
%
% Initialize chart after initial correction and add to curve segment point
% list.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

chart         = cseg.curr_chart;
chart.pt      = atlas.next_pt;
atlas.next_pt = atlas.next_pt+1;
chart.id      = chart.pt+1;
chart.dim     = atlas.dim;
chart.s       = atlas.cont.s0;
chart.bv      = atlas.cont.bv0;
chart.P.mark  = zeros(1,numel(chart.bv));
chart.remesh  = 0;

prob = save_funcs(prob, chart.id);

[prob, cseg]  = cseg.add_chart(prob, chart);
chart = cseg.ptlist{1};

% projection
chart.ics = prob.efunc.p_idx;
chart.xp  = chart.x(chart.ics);
TSp       = chart.TS(chart.ics,:);
try
  chart.TSp = cseg.orth(TSp);
catch ME
  if (strcmp(ME.identifier,'MATLAB:square'))
    msg = ['The initial chart is located at the singular point (', ...
      num2str(chart.xp'),') of the projected geometry.'];
    causeException = MException('MATLAB:init_atlas:singularity',msg);
    ME = addCause(ME,causeException);
  end
  rethrow(ME)
end
chart.G   = (TSp' * TSp) \ (TSp' * chart.TSp);

chart.rmProps   = union(chart.rmProps, {'xp' 'TSp' 'G'});
cseg.ptlist{1}  = chart;
cseg.curr_chart = cseg.ptlist{1};

end
