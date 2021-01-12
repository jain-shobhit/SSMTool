function [prob, atlas, cseg] = init_admissible(atlas, prob, cseg, S)
%INIT_ADMISSIBLE  Remove inadmissible directions.
%
% Restrict continuation to admissible directions from the initial chart.
% Flush and terminate if no admissible directions.
%
% Identical to atlas1d_v5.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: init_admissible.m 3087 2019-04-04 19:54:09Z hdankowicz $

flags = (S.dir_flags==atlas.IsAdmissible);
if S.ep_flag==1 && any(flags)
  chart            = cseg.curr_chart;
  chart.v(~flags)  = 0.5*(atlas.cont.Rmarg*chart.R); % Make vertices interior to sphere
  chart.bv(~flags) = [];                             % Make directions unavailable
  cseg.curr_chart  = chart;
else
  [prob, atlas, cseg] = init_admissible@AtlasBase(atlas, prob, cseg, S);
end

end
