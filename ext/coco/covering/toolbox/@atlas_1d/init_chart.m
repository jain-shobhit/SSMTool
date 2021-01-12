function [prob, atlas, cseg] = init_chart(atlas, prob, cseg)
% Initialize initial chart
chart      = cseg.curr_chart;
chart.pt   = 0;
chart.s    = [-1 1];
chart.im   = {'prepend' 'append'};
atlas.PtMX = abs(atlas.cont.PtMX);
if sum(atlas.PtMX)==0
  chart.s    = 1;
  chart.im   = {'append'};
  atlas.PtMX = 0;
elseif any(atlas.PtMX==0)
  idx        = atlas.PtMX~=0;
  chart.s    = chart.s(idx);
  chart.im   = chart.im(idx);
  atlas.PtMX = atlas.PtMX(idx);
end
% this update ensures smaller angles between tangents at first curve
% segment, in particular, with adaptation enabled
cseg.src_chart  = chart;
prob            = coco_emit(prob, 'update', cseg);
[prob, cseg]    = cseg.add_chart(prob, chart);
cseg.curr_chart = cseg.ptlist{1};

CurveSegment.log_angle(prob, 2, atlas.cont.LogLevel, ...
  cseg.curr_chart.t, prob.efunc.p_idx, mfilename, 'init_chart');
end
