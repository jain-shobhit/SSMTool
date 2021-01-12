function isdata = coco_is_chart_data(chart, fid)
isdata  = any(strcmp(fid, chart.private.data(:,1)));
end