function [prob, atlas, cseg] = merge(atlas, prob, cseg)
% MERGE   Merge a chart into the atlas
%
% Adds a new chart to chart tree and adjust polygonal representations

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

chart1 = cseg.ptlist{end};
edge   = atlas.edge;
charts = atlas.charts;

if isempty(atlas.tree)
  atlas.tree = create_tree(numel(chart1.xp)); % Initialize the tree
else
  % Construct list of all overlapping charts
  list = list_intersecting_charts(atlas.tree, -1, chart1.xp, chart1.R);
  chart1.intsct = list;
  
  for i=1:numel(list)
    chart2 = charts{list(i)};
    if ~atlas.isclose(chart1, chart2)
      continue
    end
    edge   = edge + 1;
    chart1 = chop_polyhedron(edge, chart1, chart2);
    if chart2.P.n~=0
      edge   = edge + 1;
      chart2 = chop_polyhedron(edge, chart2, chart1);
    end
    chart2.intsct   = [chart2.intsct chart1.id];
    charts{list(i)} = chart2; % Update existing charts
  end
  if charts{1}.ep_flag==1
    charts{1}.bv = [];
  end
end

charts = [charts chart1]; % Append new chart to atlas

border      = setdiff(union(chart1.intsct, chart1.id), 1);
brdr_charts = charts(border);
idx         = cellfun(@(x) x.ep_flag==1, brdr_charts);
border      = border(idx);
for i=1:numel(border)
  chart = charts{border(i)};
  chart.bv = [];
  charts{border(i)} = chart;
end

atlas.boundary = union([atlas.boundary; chart1.id], border);
bdry_charts    = charts(atlas.boundary);
idx            = cellfun(@(x) ~isempty(x.bv), bdry_charts);
atlas.boundary = atlas.boundary(idx);

if isempty(atlas.boundary) || (atlas.next_pt>atlas.cont.PtMX)
  chart = cseg.ptlist{end};
  chart.pt_type    = 'EP';
  chart.ep_flag    = 1;
  cseg.ptlist{end} = chart;
end

atlas.charts = charts;
atlas.edge   = edge;
atlas = atlas.add_chart2tree(chart1.id); % Insert new chart in tree

end
