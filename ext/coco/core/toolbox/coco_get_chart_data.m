function data = coco_get_chart_data(chart, fid)
idx  = strcmp(fid, chart.private.data(:,1));
data = chart.private.data{idx,2};
end
