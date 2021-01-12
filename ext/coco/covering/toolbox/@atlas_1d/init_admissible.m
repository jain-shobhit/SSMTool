function [prob, atlas, cseg] = init_admissible(atlas, prob, cseg, S)
if S.ep_flag == 1
  % compute subset of admissible directions
  flags = (S.dir_flags==atlas.IsAdmissible);
  cseg.curr_chart.s  = cseg.curr_chart.s(flags);
  cseg.curr_chart.im = cseg.curr_chart.im(flags);
  atlas.PtMX         = atlas.PtMX(flags);
  if ~any(flags)
    coco_warn(prob, 1, prob.cont.LogLevel, ...
      '%s: %s\n%s: %s\n', mfilename, ...
      'no direction admissible', ...
      'active boundary or terminal constraints were', S.msg);
  end
else
  [prob, atlas, cseg] = init_admissible@AtlasBase(atlas, prob, cseg, S);
end
end
