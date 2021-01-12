function prob = save_funcs(prob, id)

[prob, funcs] = coco_save_funcs(prob); %#ok<ASGLU>
fname = fullfile(prob.run.dir, sprintf('func_data%d', id));
run     = prob.run; %#ok<NASGU>
version = {1 'full'}; %#ok<NASGU>
saved  = false;
trials = 1;
trMX   = 3;
while ~saved
  try
    save(fname, 'version', 'run', 'funcs');
    saved = true;
  catch e
    trials = trials+1;
    if trials<=trMX
      msg = sprintf('%s: caught error while saving solution file ''%s'':', mfilename, fname);
      msg = sprintf('%s\n%s', msg, e.message);
      msg = sprintf('%s\nretrying: trial %d of %d ...', msg, trials, trMX);
      coco_warn(prob, 1, 1, '%s', msg);
      pause(1);
    else
      msg = sprintf('%s: caught error while saving solution file ''%s'':', mfilename, fname);
      msg = sprintf('%s\n%s', msg, e.message);
      msg = sprintf('%s\nsaving solution file failed, giving up ...', msg);
      coco_warn(prob, 1, 1, '%s', msg);
      saved = true;
    end
  end
end

end
