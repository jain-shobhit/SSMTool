function prob = continex_add_recording(prob)
% Add extensive recording of correction data during continuation.
% Use with care, large data files will result. Use
%   data = coco_read_solution(FID, RUN, LAB, 'data');
% to access all recorded data. Total number of correction steps (corr.its)
% and number of trials (cseg.trial) are stored in BD.

data = data_init();

data = coco_func_data(data);

prob = coco_add_slot(prob, 'continex.rec', @corr_begin,  data, 'corr_begin' );
prob = coco_add_slot(prob, 'continex.rec', @corr_step,   data, 'corr_step'  );
prob = coco_add_slot(prob, 'continex.rec', @corr_end,    data, 'corr_end'   );
prob = coco_add_slot(prob, 'continex.rec', @corr_sample, data, 'corr_sample');

prob = coco_add_slot(prob, 'continex.rec', @save_red, data, 'save_reduced');

prob = coco_add_slot(prob, 'continex.rec', @bd_add, data, 'bddat');
end

function data = data_init(data)
data.trial   = 0;
data.it      = 0;
data.history = cell(1,0);
end

function data = corr_begin(prob, data, varargin)
% varargin = { method[='nwtn'] chart x }
data.trial = data.trial+1;
data       = add_history(prob, data, 'corr.begin', varargin{2:3});
end

function [data stop msg] = corr_step(prob, data, varargin)
% varargin = { method[='nwtn'] chart x }
data.it = data.it+1;
data    = add_history(prob, data, 'corr.step', varargin{2:3});
stop    = false;
msg     = '';
end

function data = corr_end(prob, data, varargin)
% varargin = { method[='nwtn' flag[=accept] chart x }
%          | { method[='nwtn' flag[=fail|stop]
flag = varargin{2};
id   = sprintf('corr.%s', flag);
switch flag
  case 'accept'
    data = add_history(prob, data, id, varargin{3:4});
  otherwise
    data = add_history(prob, data, id, struct(), []);
end
end

function data = corr_sample(prob, data, varargin)
% varargin = { x }
data = add_history(prob, data, 'corr.sample', struct(), varargin{1});
end

function [data res] = save_red(prob, data, varargin) %#ok<INUSL>
% varargin = { chart [chart1] }
res  = data.data;
data = data_init(data);
end

function [data res] = bd_add(prob, data, cmd, varargin) %#ok<INUSL>
% varargin = [cmd='init'] { }
%          | [cmd='data'] { chart }
switch cmd
  case 'init'
    res = { 'cseg.trial' 'corr.its' };
  case 'data'
    res = { data.trial data.it };
end
end

function data = add_history(prob, data, id, chart, x)
if coco_exist('cseg.prcond', 'func', prob)
  pr_data = coco_get_func_data(prob, 'cseg.prcond', 'data');
  data.history = [
    data.history
    { id , data.trial , data.it , ...
    pr_data.data, prob.corr , prob.cseg , chart, x }
    ];
else
  data.history = [
    data.history
    { id , data.trial , data.it , ...
    struct(), prob.corr , prob.cseg , chart, x }
    ];
end
end
