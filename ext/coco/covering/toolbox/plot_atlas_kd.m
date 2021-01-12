function plot_atlas_kd(charts, varargin)
% PLOT_ATLAS_KD   Skeletal visualization of one- and two-dimensional atlases
%
% PLOT_ATLAS_KD(ATLAS, VARARGIN)
%
% ATLAS - Manifold atlas.
%
% VARARGIN = {(DIM | ICS1, ICS2, ICS3), OPT}
% OPT      = ('polyhedra' | 'basepoints' | 'spheres')
% DIM   - Embedding dimension.
% ICS1, ICS2, ICS3 - active continuation parameter indices for graphing

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

if ischar(varargin{end})
  assert(any(strcmpi(varargin{end}, ...
    {'polyhedra', 'basepoints', 'spheres'})), 'Unsupported plotting type');
  opt = varargin{end};
  varargin(end) = [];
else
  opt = 'polyhedra';
end
num=0;

dim = length(varargin);
switch dim
  case 1
    dim = varargin{1};
    ics1 = 1;
    ics2 = 2;
    ics3 = 3;
  case 2
    [ics1, ics2] = deal(varargin{:});
  case 3
    [ics1, ics2, ics3] = deal(varargin{:});
  otherwise
    error('%s: wrong dimension', mfilename);
end

hold on
switch dim
  case 2
    for i=1:numel(charts)
      chart = charts{i};
      pt = chart.x(chart.ics([ics1 ics2]));
      if strcmpi(opt, 'basepoints')
        plot(pt(1), pt(2), 'k.', 'MarkerSize', 12);
        text(pt(1), pt(2), sprintf('%d',chart.id), 'FontSize', 10);
      end
      if strcmpi(opt, 'spheres')
        pts = repmat(chart.x, [1, 2]) + chart.TS*chart.G*chart.R*[1 -1];
        pts = pts(chart.ics([ics1 ics2]),:);
        plot(pts(1,:), pts(2,:), 'ro');
      end
      if chart.P.n>0 && strcmpi(opt, 'polyhedra')
        pts = repmat(chart.x, [1, 2]) + ...
          chart.TS*chart.G*[ chart.P.v{1} chart.P.v{2} ];
        pts = pts(chart.ics([ics1 ics2]),:);
        if isempty(chart.bv)
          plot(pts(1,:), pts(2,:), 'r');
        else
          plot(pts(1,:), pts(2,:), 'b', 'LineWidth', 1.5);
        end
      end
    end
  case 3
    for i=1:numel(charts)
      chart = charts{i};
      pt = chart.x(chart.ics([ics1 ics2 ics3]));
      if strcmpi(opt, 'basepoints')
        plot3(pt(1), pt(2), pt(3), 'k.', 'MarkerSize', 6);
        text(pt(1), pt(2), pt(3), sprintf('%d',chart.id), 'FontSize', 5);
      end
      if strcmpi(opt, 'spheres')
        switch size(chart.TS,2)
          case 2
            th=0:0.1:2*pi;
            pts = repmat(chart.x, [1, numel(th)]) + ...
              chart.TS*chart.G*chart.R*[cos(th); sin(th)];
            pts = pts(chart.ics([ics1 ics2 ics3]),:);
            if isempty(chart.bv)
              patch('Faces', 1:numel(th), 'Vertices', pts', ...
                'FaceColor', 'green', 'FaceAlpha', 0.4)
            else
              patch('Faces', 1:numel(th), 'Vertices', pts', ...
                'FaceColor', 'blue', 'FaceAlpha', 0.4)
            end
            plot3(pts(1,:), pts(2,:), pts(3,:), 'k', 'LineWidth', 1);
          case 1
            pts = repmat(chart.x, [1, 2]) + ...
              chart.TS*chart.G*chart.R*[1 -1];
            pts = pts(chart.ics([ics1 ics2 ics3]),:);
            plot3(pts(1,:), pts(2,:), pts(3,:), 'k.');
          otherwise
        end
      end
      if chart.P.n>0 && strcmpi(opt, 'polyhedra')
        num=num+1;
        switch size(chart.TS,2)
          case 2
            v = cell2mat(chart.P.v);
            vrts = sortvertices(chart.P);
            pts = repmat(chart.x, [1, numel(vrts)]) + ...
              chart.TS*chart.G*v(:, vrts);
            pts = pts(chart.ics([ics1 ics2 ics3]),:);
            if isempty(chart.bv)
              patch('Faces', 1:chart.P.n, 'Vertices', pts', ...
                'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1)
            else
              patch('Faces', 1:chart.P.n, 'Vertices', pts', ...
                'FaceColor', 'blue', 'FaceAlpha', 0.4)
            end
          case 1
            pts = repmat(chart.x, [1, 2]) + ...
              chart.TS*chart.G*chart.R*[1 -1];
            pts = pts(chart.ics([ics1 ics2 ics3]),:);
            plot3(pts(1,:), pts(2,:), pts(3,:), 'k', 'LineWidth', 1);
          otherwise
        end
      end
    end
  otherwise
end
hold off
%num

end

function vs = sortvertices(P)

temp = cell2mat(P.faceV');
vs   = temp(1,:);
lookfor   = temp(1,2);
temp(1,:) = [0, 0];
for k=2:P.nFaces-1
  i = find(lookfor==temp(:,1));
  if isempty(i)
    i = find(lookfor==temp(:,2));
    lookfor = temp(i(1),1);
  else
    lookfor = temp(i(1),2);
  end
  vs(k+1) = lookfor;
  temp(i(1),:) = [0, 0];
end

end
