function [prob, atlas, cseg] = flush(atlas, prob, cseg)
% Flush chart to disk and screen output.

% flush point list
[prob, atlas, cseg] = atlas.flush@AtlasBase(prob, cseg);

if cseg.Status == cseg.CurveSegmentOK
  chart = cseg.ptlist{end};
  
  % flush last point into chart_list
  atlas.chart_list{1} = chart;
  
  % end of branch if It>=PtMX.
  if isempty(atlas.PtMX) || (chart.pt>=atlas.PtMX(1))
    cseg.Status = cseg.BoundaryPoint;
  end
end

% Handle end of branch.
if cseg.Status~=cseg.CurveSegmentOK && numel(atlas.chart_list)>1
  atlas.PtMX       = atlas.PtMX(2:end);
  atlas.chart_list = atlas.chart_list(2:end);
  chart            = atlas.chart_list{1};
  atlas.func_list  = atlas.func_list(2:end);
  atlas.first      = true;
  prob = coco_restore_funcs(prob, atlas.func_list{1});
  prob = AtlasBase.bddat_set(prob, 'ins_mode', chart.im);
  atlas.PrintHeadLine = true;
  cseg.Status=cseg.CurveSegmentOK;
end

end
