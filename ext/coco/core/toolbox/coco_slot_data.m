function varargout = coco_slot_data(fid, results)

idx  = strcmp(fid, results(:,1));
if any(idx)
  [varargout{1:nargout}] = deal(results{idx,2:nargout+1});
else
  [varargout{1:nargout}] = deal([]);
end

end
