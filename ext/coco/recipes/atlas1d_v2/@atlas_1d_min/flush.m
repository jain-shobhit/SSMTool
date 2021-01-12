function [prob atlas cseg] = flush(atlas, prob, cseg)
%FLUSH   Flush point list.
%
% Use superclass flush function to write point list to screen and disk.
% Terminate when number of successfully completed continuation steps equals
% PtMX.
%
% Identical to atlas1d_v1.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: flush.m 2839 2015-03-05 17:09:01Z fschild $

[prob atlas cseg] = atlas.flush@AtlasBase(prob, cseg);
if cseg.Status==cseg.CurveSegmentOK
  atlas.base_chart = cseg.ptlist{end};
  if atlas.base_chart.pt>=atlas.cont.PtMX
    cseg.Status = cseg.BoundaryPoint;
  end
end

end
