function opts = coco_save_reduced(opts, chart, chart1)
%COCO_SAVE_REDUCED  Save reduced amount of data for labelled solution.
%
%   OPTS = COCO_SAVE_REDUCED(OPTS) saves a reduced amount of data for a
%   labelled solution that is not a restart solution. This data should be a
%   subset of the data saved with CONT_SAVE_FULL. This function calls
%   OPTS.CONT.SAVE_REDUCED to obtain a list of additional class names that
%   should be saved. If this list is empty (default), this function saves
%   nothing.
%
%   See also: cont_save_full, coco_set, coco_opts
%

if nargin<=2
  chart1 = struct();
end

fname = fullfile(opts.run.dir, sprintf('sol%d', chart.lab));
run   = opts.run; %#ok<NASGU>
[opts fids data] = coco_emit(opts, 'save_reduced', chart, chart1);
data  = [ fids data ]; %#ok<NASGU>
version = {1 'reduced'}; %#ok<NASGU>

saved  = false;
trials = 1;
trMX   = 3;
while ~saved
  try
    save(fname, 'version', 'run', 'data');
    saved = true;
  catch e
    trials = trials+1;
    if trials<=trMX
      msg = sprintf('%s: caught error while saving solution file ''%s'':', mfilename, fname);
      msg = sprintf('%s\n%s', msg, e.message);
      msg = sprintf('%s\nretrying: trial %d of %d ...', msg, trials, trMX);
      coco_warn(opts, 1, 1, '%s', msg);
      pause(1);
    else
      msg = sprintf('%s: caught error while saving solution file ''%s'':', mfilename, fname);
      msg = sprintf('%s\n%s', msg, e.message);
      msg = sprintf('%s\nsaving solution file failed, giving up ...', msg);
      coco_warn(opts, 1, 1, '%s', msg);
      saved = true;
    end
  end
end
