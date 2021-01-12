function [prob atlas cseg correct] = init_prcond(atlas, prob, chart)
%INIT_PRCOND   Initialize projection condition.
%
% Initialize curve segment projection condition and proceed to correct if
% chart.t is empty.
%
% Identical to atlas1d_v1.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_prcond.m 2839 2015-03-05 17:09:01Z fschild $

chart.R       = 0;
chart.pt      = -1;
chart.pt_type = 'IP';
chart.ep_flag = 1;
[prob cseg]   = CurveSegment.create_initial(prob, chart);
correct       = cseg.correct;

end
