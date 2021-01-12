function [prob, atlas, cseg] = init_admissible(atlas, prob, cseg, S)
%INIT_ADMISSIBLE  Remove inadmissible directions.
%
% Restrict continuation to admissible directions from the initial chart.
% Flush and terminate if no admissible directions.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

flags = (S.dir_flags==atlas.IsAdmissible);
if S.ep_flag==1 && any(flags)
  chart = cseg.curr_chart;
  chart.P.mark(~flags) = 1;
  cseg.curr_chart  = chart;
else
  cseg.curr_chart.ep_flag = 0;
  [prob, atlas, cseg] = init_admissible@AtlasBase(atlas, prob, cseg, S);
end

end
