function [prob, atlas, cseg, correct] = init_prcond(atlas, prob, chart)
% Initialize initial chart
chart.R       = 0;
chart.pt      = -1;
chart.pt_type = 'IP';
chart.ep_flag = 1;
[prob, cseg]  = CurveSegment.create_initial(prob, chart, ...
  atlas.cont.Valpha, atlas.cont.h0, atlas.cont.NullItMX);
correct       = cseg.correct;
end
