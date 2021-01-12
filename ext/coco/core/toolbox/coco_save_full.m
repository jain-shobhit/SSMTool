function prob = coco_save_full(prob, chart, chart1)
%COCO_SAVE_FULL  Save restart data for labelled solution.
%
%   PROB = COCO_SAVE_FULL(PROB, CHART, CHART1) saves restart data for a
%   labelled solution. This data usually contains the solution, tangent
%   vector and a toolbox identifyer. This function calls
%   OPTS.CONT.SAVE_FULL to obtain a list of additional class names that
%   should be saved.
%
%   See also: coco_set, coco_opts
%

if nargin<=2
  chart1 = struct();
end

fname   = fullfile(prob.run.dir, sprintf('sol%d', chart.lab));
run     = prob.run; %#ok<NASGU>
[prob, fids1, data1] = coco_emit(prob, 'save_reduced', chart, chart1);
[prob, fids2, data2] = coco_emit(prob, 'save_full',    chart, chart1);
data    = [ fids1 data1 ; fids2 data2 ]; %#ok<NASGU>
fdata   = [ prob.efunc.identifyers' { prob.efunc.funcs(:).x_idx }' ]; %#ok<NASGU>
if isfield(prob, 'adjoint') % Harry added
  adata = [ prob.adjoint.identifyers' { prob.adjoint.funcs(:).af_idx }' ]; %#ok<NASGU>
else
  adata = cell(1,2); %#ok<NASGU>
end
version = {1 'full'}; %#ok<NASGU>

saved  = false;
trials = 1;
trMX   = 3;
while ~saved
  try
    save(fname, 'version', 'run', 'data', 'chart', 'chart1', 'fdata', 'adata');
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
