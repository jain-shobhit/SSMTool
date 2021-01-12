function [atlas prob cseg] = merge(atlas, prob, cseg)
%MERGE   Merge curve segment into atlas
%
% Check if last chart in point list is neighbor with any atlas boundary
% chart and remove overlapping directions of continuation. Reconstitute
% atlas boundary.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: merge.m 2839 2015-03-05 17:09:01Z fschild $

chart = cseg.ptlist{end};
R     = atlas.cont.h;
h     = atlas.cont.Rmarg*R;
nb    = cell(numel(chart.s),4); % Not all directions are allowed
for k=1:numel(chart.s)
  sk      = chart.s(k);
  xk      = chart.x+h*(chart.TS*sk);
  nb(k,:) = {chart, xk, sk, h};
end
for i=size(atlas.boundary,1):-1:1
  chart2 = atlas.boundary{i,1};
  if atlas.isneighbor(chart, chart2)
    x2 = atlas.boundary{i,2};
    if norm(chart.TS'*(x2-chart.x))<R
      atlas.boundary(i,:) = [];
    end
    for k=size(nb,1):-1:1
      x1 = nb{k,2};
      if norm(chart2.TS'*(x1-chart2.x))<R
        nb(k,:) = [];
      end
    end
  end
end
atlas.boundary = [nb; atlas.boundary];
for j=2:numel(cseg.ptlist)-1 % Check intermediate charts in point list
  cdata = coco_get_chart_data(cseg.ptlist{j}, atlas.cdid); % Extract chart data
  if ~isempty(cdata) && isfield(cdata, 'TS') % If branch point
    chart   = cseg.ptlist{j};
    chart.t = cdata.TS; % New tangent matrix
    [prob cseg2] = CurveSegment.create_initial(prob, chart); % Create projection condition
    chart    = cseg2.curr_chart;
    chart.TS = cseg2.prcond.TS;
    chart.pt = 0;
    chart.ignore_evs = cdata.ign; % Don't detect initial branch point
    nb = [{chart, chart.x+h*chart.TS, 1, h}; ...
      {chart, chart.x-h*chart.TS, -1, h}];
    for i=size(atlas.boundary,1):-1:1 % Merge branch switched chart into atlas
      chart2 = atlas.boundary{i,1};
      if atlas.isneighbor(chart, chart2)
        x2 = atlas.boundary{i,2};
        if norm(chart.TS'*(x2-chart.x))<R
          atlas.boundary(i,:) = [];
        end
        for k=size(nb,1):-1:1
          x1 = nb{k,2};
          if norm(chart2.TS'*(x1-chart2.x))<R
            nb(k,:) = [];
          end
        end
      end
    end
    atlas.boundary = [nb; atlas.boundary];
  end
end
if isempty(atlas.boundary)
  chart.pt_type    = 'EP';
  chart.ep_flag    = 1;
  cseg.ptlist{end} = chart;
end

end
