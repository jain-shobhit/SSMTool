function atlas = update_atlas(atlas, ID)
% UPDATE_ATLAS   Update the polyhedral mesh
%
% ATLAS = UPDATE_ATLAS(ATLAS, ID)
%
% ATLAS - Manifold atlas
% ID    - Integer id of chart with new radius

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

edge   = atlas.edge;
charts = atlas.charts;

chart1 = charts{ID};
nOut   = chart1.intsct;
nIn    = [];
while ~isempty(nOut)
  chart2 = charts{nOut(1)};
  chart2 = create_polyhedron(chart2, atlas.cont.Rmarg);

  % only merge with merged charts
  chart2.intsct = setdiff(chart2.intsct, ID);
  intsct = setdiff(chart2.intsct, nOut);
  for j=1:numel(intsct)
    chart3 = charts{intsct(j)};
    if ~atlas.isclose(chart2, chart3)
      continue
    end
    edge   = edge + 1;
    chart2 = chop_polyhedron(edge, chart2, chart3);
    if ismember(chart3.id, nIn) && chart3.P.n~=0
      edge   = edge + 1;
      chart3 = chop_polyhedron(edge, chart3, chart2);
    end
    charts{intsct(j)} = chart3;
  end
  charts{nOut(1)} = chart2;
  nIn = [nIn nOut(1)]; %#ok<AGROW>
  nOut(1) = [];
end

% re-insert resized chart
chart1 = create_polyhedron(chart1, atlas.cont.Rmarg);

list   = list_intersecting_charts(atlas.tree, ID, chart1.xp, chart1.R); % Construct list of all overlapping charts
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
  charts{list(i)} = chart2;
end

charts{ID} = chart1;

border      = union([nIn ID], chart1.intsct);
brdr_charts = charts(border);
idx         = cellfun(@(x) x.ep_flag==1, brdr_charts);
border      = border(idx);
for i=1:numel(border)
  chart = charts{border(i)};
  chart.bv = [];
  charts{border(i)} = chart;
end

atlas.boundary = union(atlas.boundary, [nIn ID]');
bdry_charts    = charts(atlas.boundary);
idx            = cellfun(@(x) ~isempty(x.bv), bdry_charts);
atlas.boundary = atlas.boundary(idx);
atlas.charts   = charts;
atlas.edge     = edge;

end
