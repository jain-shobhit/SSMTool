function coco_plot_u(varargin)
if ishandle(varargin{1})
  [h opts x fmt] = coco_deal(varargin{:}, '.-');
else
  [h opts x fmt] = coco_deal(gca, varargin{:}, '.-');
end

list = coco_get_func_data(opts, 'efunc', 'x0idx');
xticks = [];
xlabs  = {};
for i=1:size(list,1)
  xlabs  = [ xlabs , list{i,1} ]; %#ok<AGROW>
  xticks = [ xticks, list{i,2}(1) ]; %#ok<AGROW>
end
plot(h, x, fmt);
set(h, 'XTick', xticks, 'XTickLabel', xlabs, 'XLim', [0 numel(x)+1]);
grid(h, 'on');
end
