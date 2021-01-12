function coco_print_workspace
% COCO_PRINT_WORKSPACE Print all workspace variables of caller to screen.

dstat = get(0, 'Diary');
dfile = get(0, 'DiaryFile');

try
  diary off
  if exist('dump.txt', 'file')
    delete('dump.txt');
  end
  diary('dump.txt');
  diary on
  fprintf('list of workspace variables, see ''help whos'' for format\n');
  s = dbstack;
  if numel(s)<=1
    fprintf('workspace is ''base''\n');
  else
    fprintf('workspace of function %s at line %d\n', s(2).name, s(2).line);
  end
  s = evalin('caller', 'whos()');
  fprintf('%d variable(s) defined\n\n', numel(s));
  for i=1:numel(s)
    fprintf('**********************************************************************\n');
    fprintf('variable %d [%s] = \n\n', i, s(i).name);
    disp(s(i));
    evalin('caller', sprintf('display(%s)', s(i).name));
  end
  fprintf('**********************************************************************\n');
catch e
  msg = sprintf('%s: cought error while writing dump file:', mfilename);
  msg = sprintf('%s\n%s', msg, e.message);
  disp(msg);
end

diary off
set(0, 'DiaryFile', dfile);
set(0, 'Diary', dstat);
end
