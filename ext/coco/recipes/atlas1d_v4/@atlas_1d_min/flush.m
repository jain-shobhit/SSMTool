function [prob atlas cseg] = flush(atlas, prob, cseg)
%FLUSH   Flush point list.
%
% Merge last successfully located chart into atlas. Use superclass flush
% function to write point list to screen and disk. Terminate when number of
% successfully completed continuation steps equals PtMX along both
% directions of continuation or when the atlas boundary is empty.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: flush.m 2839 2015-03-05 17:09:01Z fschild $

if cseg.Status==cseg.CurveSegmentOK
  [atlas cseg] = atlas.merge(cseg);
end
[prob atlas cseg] = atlas.flush@AtlasBase(prob, cseg);
closed = isempty(atlas.boundary);
if cseg.Status==cseg.CurveSegmentOK
  if closed || (atlas.boundary{1,1}.pt>=atlas.cont.PtMX)
    cseg.Status = cseg.BoundaryPoint;
  end
elseif cseg.Status==cseg.BoundaryPoint && ~closed
  atlas.boundary{1,1}.pt=atlas.cont.PtMX;
  atlas.boundary = atlas.boundary([2:end 1],:);
  prob = AtlasBase.bddat_set(prob, 'ins_mode', 'append');
  if atlas.boundary{1,1}.pt<atlas.cont.PtMX
    cseg.Status = cseg.CurveSegmentOK;
  end
end

end
