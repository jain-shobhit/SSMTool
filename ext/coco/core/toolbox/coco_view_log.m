function coco_view_log(runid, mode)

if nargin<2
  if usejava('desktop')
    mode = 'edit';
  else
    mode = 'type';
  end
end

logfname = coco_fname(runid, 'coco_log.txt');
switch mode
  case 'edit'
    edit(logfname);
  otherwise
    type(logfname);
end

end
