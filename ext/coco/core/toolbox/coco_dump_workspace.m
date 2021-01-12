function coco_dump_workspace
% COCO_DUMP_WORKSPACE Save all workspace variables of caller to dump file.

% bug: check for version and save additionally in format that can be loaded in R2009b
try
  evalin('caller', 'save dump');
catch e
  msg = sprintf('%s: cought error while writing dump file:', mfilename);
  msg = sprintf('%s\n%s', msg, e.message);
  disp(msg);
end
end
