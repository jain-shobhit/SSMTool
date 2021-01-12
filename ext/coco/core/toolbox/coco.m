function bd = coco(varargin)
%COCO   Execute a continuation run.
% 
% BD = COCO(PROB, RUN, TB_SEL, CONT_ARGS)
% Execute a continuation run. The resulting data will be saved under run
% name RUN. The zero problem is defined by the toolbox selected with the
% arguments passed as TB_SEL. The arguments provided as CONT_ARGS are
% passed on to the constructor of the covering toolbox.
%
% TB_SEL = { TB_NAME ISOL_T SOL_T TB_ARGS }
% Use zero problem from toolbox TB_NAME. The constructor for the zero
% problem is the function with name sprintf('%s_%s2%s', TB_NAME, ISOL_T,
% SOL_T) and constructs a zero problem starting at a solution of type
% ISOL_T and computing a branch of solutions of type SOL_T. Arguments
% required by the problem constructor are passed as TB_ARGS.
%
% TB_SEL = { @TB_CTOR TB_ARGS }
% Use zero problem constructed by the problem constructor TB_CTOR. Arguments
% required by the problem constructor are passed as TB_ARGS.
%
% TB_SEL = { [] }
% Use zero problem defined in PROB.
%
% CONT_ARGS = { [DIM=1] [PARS={} [PINT={}]] }
% The syntax for CONT_ARGS is defined by the constructor of the covering
% algorithm. The default is covering_create, which accepts up to 3 input
% arguments. DIM is the dimension of the solution manifold, PARS a cell
% array of strings with the names of parameters and PINT a cell array of
% intervals defining computational boundaries for the parameters listed in
% PARS. An empty interval may be used to indicate the absence of
% computational boundaries.
%
% See also: COVERING_CREATE, COCO_BD_READ, COCO_BD_PRINT, COCO_PROB,
% COCO_SET

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco.m 2974 2017-01-10 21:14:29Z hdankowicz $

%% clean up after interrupted or crashed session

cleanup = coco_cleanup();
cleanup.call(@coco_clear_cache);
ptrs = coco_func_data.pointers('copy');
cleanup.call(@coco_func_data.pointers, 'set', ptrs);

%% process input arguments
%  varargin = { [prob], runid, [], ... }
%  varargin = { [prob], runid, @tbxctor, ... }
%  varargin = { [prob], runid, toolbox, fromST, toST, ... }

p = coco_stream(varargin{:});

% bug: make argument prob compulsory
if isstruct(p.peek)
	prob = p.get;
elseif isempty(p.peek)
  p.skip;
  prob = coco_set();
else
  prob = coco_set();
end

runid   = p.get;
Toolbox = p.get;

if isempty(Toolbox)
  Toolbox = 'empty';
  TbxCtor = @empty_ctor;
  from_st = '';
  to_st   = '';
elseif isa(Toolbox, 'function_handle')
  TbxCtor = Toolbox;
  Toolbox = func2str(Toolbox);
  from_st = '';
  to_st   = '';
else
  from_st = p.get;
  to_st   = p.get;
  cfname  = sprintf('%s_%s2%s', Toolbox, from_st, to_st);
  TbxCtor = str2func(cfname);
end

%% create and clean data directory
opts = get_settings(prob);

if ischar(runid)
  data_dir = fullfile(opts.data_dir, runid);
else
  data_dir = fullfile(opts.data_dir, runid{:});
end

if ~exist(data_dir, 'dir')
	[status,msg,msgid]=mkdir(data_dir); %#ok<NASGU,ASGLU>
  if ~exist(data_dir, 'dir')
    error('%s: could not create directory ''%s''\n%s: %s', ...
      mfilename, data_dir, mfilename, msg);
  end
end

if opts.CleanData
	delete(fullfile(data_dir, '*.mat'));
	delete(fullfile(data_dir, 'coco_*.txt'));
end

%% create run data and message logging information

run.runid     = runid;
run.toolbox   = Toolbox;
run.tbxctor   = TbxCtor;
run.isol_type = from_st;
run.sol_type  = to_st;
run.dir       = data_dir;
run.bdfname   = fullfile(data_dir, 'bd.mat');
run.logname   = fullfile(data_dir, 'coco_log.txt');
run.loghan    = fopen(run.logname, 'w');
cleanup.fclose(run.loghan);
run.scrname   = fullfile(data_dir, 'coco_scr.txt');
run.scrhan    = fopen(run.scrname, 'w');
cleanup.fclose(run.scrhan);
if opts.LogLevel(1)>0
  run.logPrioMN = opts.LogLevel(1);
  if numel(opts.LogLevel)>=2
    run.scrPrioMN = opts.LogLevel(2);
  else
    run.scrPrioMN = 1;
  end
else
  run.logPrioMN = 0;
  run.scrPrioMN = 0;
end
prob.run = run;

coco_log(prob, 2, opts.LogLevel, ...
  '%s: entering ''coco'', start building problem\n', mfilename);

coco_log(prob, 1, opts.LogLevel, ...
  'MATLAB version %s on architecture %s\n\n', version, computer('arch')) 

coco_log(prob, 1, opts.LogLevel, 'run: { ');
run_fields = {'runid' 'toolbox' 'tbxctor' 'isol_type' 'sol_type' ...
  'dir' 'bdfname' 'logname' 'scrname'};
for i=1:numel(run_fields)
  run_field = run_fields{i};
  val = run.(run_field);
  if ischar(val)
    coco_log(prob, 1, opts.LogLevel, '%s=''%s'' ', run_field, val);
  elseif isinteger(val)
    coco_log(prob, 1, opts.LogLevel, '%s=%d ', run_field, val);
  elseif isa(val, 'function_handle')
    coco_log(prob, 1, opts.LogLevel, '%s=@%s ', run_field, func2str(val));
  elseif iscellstr(val)
    coco_log(prob, 1, opts.LogLevel, '%s={ ', run_field);
    for k=1:numel(val)
      coco_log(prob, 1, opts.LogLevel, '''%s'' ', val{k});
    end
    coco_log(prob, 1, opts.LogLevel, '} ');
  end
end
coco_log(prob, 1, opts.LogLevel, '}\n\n');

prob = coco_add_signal(prob, 'save_bd', 'coco');
prob = coco_add_slot(prob, 'run', @save_run, [], 'save_bd');

%% call constructor of toolbox
coco_log(prob, 2, opts.LogLevel, ...
  '%s: calling constructor of top-level toolbox ''%s''\n', mfilename, Toolbox);
prob = run.tbxctor(prob, '', p);

%% call constructor of adjoint
if isfield(prob, 'adjoint')
  prob = coco_add_func_after(prob, 'efunc', @adjoint_add);
end

%% call constructor of continuer
coco_log(prob, 2, opts.LogLevel, ...
  '%s: calling constructor of covering toolbox ''%s''\n', mfilename, func2str(opts.ContAlg));
prob = opts.ContAlg(prob, p);
prob = coco_add_slot(prob, 'bd', @save_bd, [], 'save_bd');
prob = coco_add_slot(prob, 'bddat', @save_bddat, [], 'save_bd');

%% check that all arguments were used
if numel(p)>0
	error('%s: too many arguments', mfilename);
end

%% call entry function of continuer
%  the call to save is necessary to eliminate side-effects caused by
%  using handle-classes like coco_func_data as part of the prob structure
coco_log(prob, 2, opts.LogLevel, ...
  '%s: construction finished\n\n', mfilename);
fhan = coco_log(prob, 1, opts.LogLevel);
if ~isempty(fhan)
  coco_print_opts(fhan, prob);
  fprintf(fhan, '\n');
  coco_print_funcs(fhan, prob);
  fprintf(fhan, '\n');
  coco_print_slots(fhan, prob);
  fprintf(fhan, '\n');
  coco_print_sigs(fhan, prob);
  fprintf(fhan, '\n');
end
coco_log(prob, 2, opts.LogLevel, ...
  '%s: entering finite state machine\n', mfilename);
coco_log(prob, 1, opts.LogLevel, '%s\n', ...
  '*********************************************************************');
try
  % prob = coco_save_funcs(prob);
  prob = prob.cont.run(prob);
  coco_log(prob, 2, opts.LogLevel, ...
    '\n%s: computation finished successfully\n', mfilename);
catch err
  coco_log(prob, 1, opts.LogLevel, ...
    '\n%s: computation finished with error\n', mfilename);
  coco_log(prob, 1, opts.LogLevel, 'error: %s\n', err.message);
  coco_log(prob, 1, opts.LogLevel, 'stack:\n');
  stack = err.stack;
  for i=1:numel(stack)
    coco_log(prob, 1, opts.LogLevel, '  %s:%d:<%s>\n', ...
      stack(i).name, stack(i).line, stack(i).file);
  end
  rethrow(err);
end

%% return bifurcation diagram
if nargout==1
  bd = prob.bd;
end

end

%%
function opts = get_settings(prob)
opts             = coco_get(prob, 'coco');
defaults.ContAlg = 'covering';
opts             = coco_merge(defaults, opts);
if ischar(opts.ContAlg)
  cfname         = sprintf('%s_create', opts.ContAlg);
  opts.ContAlg   = str2func(cfname);
end
end

%%
function prob = empty_ctor(prob, fid, varargin) %#ok<INUSD>
end

function [data res] = save_run(prob, data)
res = prob.run;
end

function [data res] = save_bd(prob, data)
res = prob.bd;
end

function [data res] = save_bddat(prob, data)
res = prob.bddat;
end
