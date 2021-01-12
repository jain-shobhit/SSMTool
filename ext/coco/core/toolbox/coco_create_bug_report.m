function coco_create_bug_report(scriptname)
if ~(scriptname(end-1)=='.' && scriptname(end)=='m')
  scriptname = [scriptname '.m'];
end
scriptpath = fullfile(pwd, scriptname);
if exist(scriptpath, 'file')==2
  dumpfile = fullfile(pwd, 'dump.txt');
  fhan = fopen(dumpfile, 'w');
  try
    msg = sprintf('%s: running script ''%s''', mfilename, scriptname);
    msg = sprintf('%s\nfull path: <%s>', msg, scriptpath);
    print(fhan, msg);
    run(scriptpath);
    msg = sprintf('%s: no error occurred', mfilename);
    print(fhan, msg);
    fclose(fhan);
  catch e
    tos = e.stack(1);
    msg = sprintf('%s: cought error ''%s'' in function ''%s'' at\nline %d of file <%s>', ...
      mfilename, e.identifier, tos.name, tos.line, tos.file);
    msg = sprintf('%s\nerror message: %s\nrerunning script to record information for bug report', ...
      msg, e.message);
    print(fhan, msg);
    fclose(fhan);
    
    if strcmpi(tos.file,scriptpath)
      % bug: check that error occurs in coco package script, see comments below
      new_scriptpath = fullfile(pwd, sprintf('dump_%s', scriptname));
      check_file(new_scriptpath);
      src = fopen(scriptpath, 'r');
      dst = fopen(new_scriptpath, 'w');
      for i=1:tos.line-1
        fprintf(dst, '%s\n', fgetl(src));
      end
      fprintf(dst, 'coco_dump_workspace\n');
      while ~feof(src)
        fprintf(dst, '%s\n', fgetl(src));
      end
      fclose(dst);
      fclose(src);
      record_error(dumpfile, new_scriptpath);
      delete(new_scriptpath);
    else
      % bug: check that error occurs in coco package function, find top
      %      level coco function where bug occurs
      % bug: enable this function for toolboxes under development (function
      %      files in current directory)
      [fpath new_file ext] = fileparts(tos.file); %#ok<ASGLU>
      new_file = fullfile(pwd, [new_file ext]);
      check_file(new_file);
      src = fopen(tos.file, 'r');
      dst = fopen(new_file, 'w');
      for i=1:tos.line-1
        fprintf(dst, '%s\n', fgetl(src));
      end
      fprintf(dst, 'coco_dump_workspace\n');
      while ~feof(src)
        fprintf(dst, '%s\n', fgetl(src));
      end
      fclose(dst);
      fclose(src);
      record_error(dumpfile, scriptpath);
      delete(new_file);
      
      fhan = fopen(dumpfile, 'a');
      msg = sprintf('%s: error information recorded', mfilename);
      msg = sprintf('%s\nPlease attach the file ''dump.mat''', msg);
      msg = sprintf('%s\nas well as *all* relevant log files ''coco_log.txt''', msg);
      msg = sprintf('%s\nto the ticket for this bug report.', msg);
      print(fhan, msg);
      fclose(fhan);
    end
  end
else
  error('%s: m-file ''%s'' not found in current directory ''%s''', ...
    mfilename, scriptname, pwd);
end
end

function check_file(fname)
assert(exist(fname, 'file')~=2, ...
  '%s: file ''%s'' exists\ndelete or rename the file and restart', mfilename, fname)
end

function print(fhan, msg)
fprintf(fhan, '%s\n', msg);
fprintf('%s\n', msg);
end

function record_error(dumpfile, scriptpath)
diary(dumpfile)
diary on
try
  run(scriptpath);
catch %#ok<CTCH>
end
diary off
end
