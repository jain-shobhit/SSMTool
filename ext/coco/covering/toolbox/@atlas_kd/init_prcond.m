function [prob, atlas, cseg, correct] = init_prcond(atlas, prob, chart)
%INIT_PRCOND   Initialize projection condition.
%
% Initialize curve segment projection condition and proceed to correct if
% chart.t is empty.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

assert(numel(prob.efunc.p_idx)>=atlas.dim, '%s: %s %d %s', mfilename, ...
  'fewer than ', atlas.dim, 'active parameters');
chart.R       = 0;
chart.pt      = -1;
chart.pt_type = 'IP';
chart.ep_flag = 1;
[prob, cseg]  = CurveSegment.create_initial(prob, chart);
correct       = cseg.correct;

end
