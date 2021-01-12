function [prob, atlas, cseg, flush] = add_chart(atlas, prob, cseg)
%ADD_CHART   Add chart to point list.
%
% Add successully located chart to curve segment point list and increment
% counter. Designate end point by 'EP' point type and flush. Terminate if
% successive charts do not qualify as "neighbors".

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

if atlas.first
  chart = cseg.ptlist{1};
  cseg.curr_chart.ignore_evs = chart.ignore_evs;
  cseg.curr_chart.pt_type    = chart.pt_type;
  cseg.curr_chart.ep_flag    = chart.ep_flag;
  cseg.curr_chart.ics        = chart.ics;
  cseg.curr_chart.remesh     = 0;
  
  prob  = save_funcs(prob, chart.id);
  [prob, cseg] = cseg.add_chart(prob, cseg.curr_chart);
  chart = cseg.ptlist{end};
  % projection
  chart.xp = chart.x(chart.ics);
  TSp      = chart.TS(chart.ics,:);
  try
    chart.TSp = cseg.orth(TSp);
  catch ME
    if (strcmp(ME.identifier,'MATLAB:square'))
      msg = ['The new chart is located at the singular point (', ...
        num2str(chart.xp'),') of the projected geometry.'];
      causeException = MException('MATLAB:init_atlas:singularity',msg);
      ME = addCause(ME,causeException);
    end
    rethrow(ME)
  end
  chart.G = (TSp' * TSp) \ (TSp' * chart.TSp);
  chart.rmProps = union(chart.rmProps, {'xp' 'TSp' 'G'});
  atlas.charts{atlas.boundary(1)} = chart;
  atlas = atlas.update_atlas(atlas.boundary(1));
  flush = false;
  return
end
chart         = cseg.curr_chart;
chart.pt      = atlas.next_pt;
atlas.next_pt = atlas.next_pt+1; % Increment chart counter
chart.id      = chart.pt+1;
chart.dim     = atlas.dim;
chart.remesh  = chart.remesh+1;

% [prob, atlas.funcs{chart.id}] = coco_save_funcs(prob);
prob = save_funcs(prob, chart.id);

[prob, cseg] = cseg.add_chart(prob, chart);

chart = cseg.ptlist{end};
% projection
chart.xp  = chart.x(chart.ics);
TSp       = chart.TS(chart.ics,:);
try
  chart.TSp = cseg.orth(TSp);
catch ME
  if (strcmp(ME.identifier,'MATLAB:square'))
    msg = ['The new chart is located at the singular point (', ...
      num2str(chart.xp'),') of the projected geometry.'];
    causeException = MException('MATLAB:init_atlas:singularity',msg);
    ME = addCause(ME,causeException);
  end
  rethrow(ME)
end
chart.G = (TSp' * TSp) \ (TSp' * chart.TSp);
chart.rmProps = union(chart.rmProps, {'xp' 'TSp' 'G'});
cseg.ptlist{end} = chart;

flush        = true;

%% step size control: update step size

chart1   = cseg.ptlist{1};
chart2   = cseg.ptlist{2};
arc_beta = subspace(chart1.TSp, chart2.TSp);

cont     = atlas.cont;
h        = chart1.h;

% This if-statement takes care of the case arc_beta==0, which produced an
% error in some versions of Matlab.
if cont.R_fac_max^2*arc_beta < cont.almax
  h_fac = cont.R_fac_max;
else
  h_fac = max(cont.R_fac_min, min(cont.R_fac_max, ...
    cont.almax/(sqrt(cont.R_fac_max)*arc_beta)));
end

h = max( min(cont.ga*h_fac*h, cont.R_max), cont.R_min);

% check residuum condition for vertices of new hypercube
[prob, chart, res] = chart_hc_res(prob, chart2, h);
while (res>cont.MaxRes) && (h>cont.R_min)
  h   = max(cont.R_min, cont.R_fac_min*h);
  [prob, chart, res] = chart_hc_res(prob, chart2, h);
end

%% Assign radius and build cube

R  = chart1.h;
ca = cos(cont.almax);
tp = chart1.t(chart.ics);
tp = tp/norm(tp);
d  = ca * (tp' * (chart2.xp-chart1.xp));
h  = max( min([h cont.R_max sqrt(R*R+d*d)]), cont.R_min);
if h<sqrt(R*R-d*d)
  chart = atlas.charts{chart1.id};
  chart.h = max(cont.R_fac_min*chart1.h, cont.R_min);
  atlas.charts{chart1.id} = chart;
  atlas.first = true;
  atlas.next_pt = atlas.next_pt - 1;
  flush    = false;
elseif ~atlas.isneighbor(chart1,chart2)
  if  chart1.h<=cont.R_min
    cseg.ptlist{end}.pt_type = 'GAP';
    cseg.ptlist{end}.ep_flag = 2;
    cseg.Status              = cseg.CurveSegmentCorrupted;
  else
    chart = atlas.charts{chart1.id};
    chart.h = max(cont.R_min, cont.R_fac_min*chart1.h);
    atlas.charts{chart1.id} = chart;
    atlas.first = true;
    atlas.next_pt = atlas.next_pt - 1;
    flush    = false;
  end
else
  if h < chart1.R
    chart = atlas.charts{chart1.id};
    chart.R = h;
    atlas.charts{chart1.id} = chart;
    atlas.tree.root = change_radius(atlas.tree.root, numel(chart.xp), ...
      chart, 1, h, false);
    atlas   = atlas.update_atlas(chart.id);
  end
  chart2.R = h;
  chart2.h = h;
  chart2   = create_polyhedron(chart2, cont.Rmarg);
  chart2.loc     = [];
  chart2.intsct  = [];
  cseg.ptlist{2} = chart2;
end

end

function [prob, chart, res] = chart_hc_res(prob, chart, h)

% Ignore residuum of arc-length constraint
prob = coco_emit(prob, 'set_mode', 'check_res');

% Compute residuum on boundary box
xx  = chart.TS*chart.G*chart.s;
res = 0;
for i=1:size(xx,2)
  [prob, chart, f] = prob.efunc.F(prob, chart, chart.x+h*xx(:,i));
  res = max(res, norm(f));
end

end
