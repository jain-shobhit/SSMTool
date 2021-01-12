function flag = coco_check_parnames(pnames)
if ~iscell(pnames)
  pnames = { pnames };
end
sargs = [ pnames(:) cell(numel(pnames),1) ]';

if nargout==0
  struct(sargs{:});
else
  flag = true;
  try
    struct(sargs{:});
  catch %#ok<CTCH>
    flag = false;
  end
end
end
