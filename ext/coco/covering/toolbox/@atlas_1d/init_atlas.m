function [prob, atlas, cseg, flush] = init_atlas(atlas, prob, cseg)
% Insert chart in point list and flush.
chart          = cseg.curr_chart;
chart.R        = atlas.cont.h0;
chart.pt_type  = 'EP';
chart.ep_flag  = 1;
cseg.ptlist{1} = chart;

% split chart into chart list
for i=1:numel(chart.s)
  chart1              = chart;
  chart1.s            = chart.s(i);
  chart1.im           = chart.im{i};
  chart1.t            = chart.s(i)*chart.t;
  [prob, chart1]      = cseg.update_p(prob, chart1);
  chart1.ignore_at    = (i>1);
  atlas.chart_list{i} = chart1;
  if i==1
    cseg.ptlist{1} = atlas.chart_list{1};
    prob = AtlasBase.bddat_set(prob, 'ins_mode', atlas.chart_list{1}.im);
  else
    [prob, atlas.func_list{i}] = coco_save_funcs(prob);
  end
end

% Flush curve segment.
flush = true;
end
