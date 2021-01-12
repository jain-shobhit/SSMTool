function [fpath flag] = coco_fname(runid, fname)

flag = true;

if ischar(runid)
  run_dir = runid;
else
  run_dir = fullfile(runid{:});
end

fpath = fullfile(run_dir, fname);
if exist(fpath, 'file')
  return
end

run_dir = fullfile('data', run_dir);
fpath   = fullfile(run_dir, fname);
if exist(fpath, 'file')
  return
end

flag = false;

if nargout<2
  if ischar(runid)
    runid_str = sprintf('''%s''', runid);
  else
    if numel(runid)>1
      runid_str = sprintf(',''%s''', runid{2:end});
    else
      runid_str = '';
    end
    runid_str = sprintf('{''%s''%s}', runid{1}, runid_str);
  end
  error('%s: file ''%s'' or run %s not found', mfilename, fname, runid_str);
end

end
