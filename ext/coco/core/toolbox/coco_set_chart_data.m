function chart = coco_set_chart_data(chart, fid, data)
idx  = strcmp(fid, chart.private.data(:,1));
chart.private.data{idx,2} = data;
end
