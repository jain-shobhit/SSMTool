function [prob, atlas, cseg] = flush(atlas, prob, cseg)
%FLUSH   Flush point list.
%
% Merge last successfully located chart into atlas. Use superclass flush
% function to write point list to screen and disk. Terminate when number of
% successfully completed continuation steps equals PtMX or when the atlas
% boundary is empty.

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

if ~isempty(cseg.ptlist)
  chart = cseg.ptlist{end};
  if ~isfield(chart, 'xp')
    prob = save_funcs(prob, chart.id);
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
    cseg.ptlist{end} = chart;
  end
end

if cseg.Status==cseg.CurveSegmentOK
  [prob, atlas, cseg] = atlas.merge(prob, cseg);
end
[prob, atlas, cseg] = atlas.flush@AtlasBase(prob, cseg);
if any(cseg.Status==[cseg.CurveSegmentOK cseg.BoundaryPoint])
  if isempty(atlas.boundary) || (atlas.next_pt>atlas.cont.PtMX)
    cseg.Status = cseg.BoundaryPoint;
  else
    cseg.Status = cseg.CurveSegmentOK;
  end
end
if cseg.Status~=cseg.CurveSegmentOK
  delete(fullfile(prob.run.dir, 'func_data*.mat'));
end

end
