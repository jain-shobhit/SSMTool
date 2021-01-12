function opts = covering_create(opts, varargin)

%% process input arguments
%  varargin = { [dim=1,] [pars={} [, pint={}]] }
p = coco_stream(varargin{:});

if numel(p)<=0
  dim   = 1;
  pars  = {};
  pints = {};
else
  if ischar(p.peek) || iscell(p.peek)
    dim = 1;
  else
    dim = p.get;
  end
  
  if numel(p)>=1
    pars = p.get;
  else
    pars = {};
  end
  
  if numel(p)>=1
    pints  = p.get;
  else
    pints  = {};
  end
end

%% get toolbox settings
[opts, cont] = get_settings(opts, dim);

%% add boundary events

if ischar(pars)
	pars  = { pars  };
end
if ~iscell(pints)
  pints = { pints };
end

for i=1:numel(pints)
  par  = pars {i};
  pint = pints{i};
  if ~isempty(pint)
    opts = coco_add_event(opts, 'EP', 'boundary', par, '<', pint(1));
    opts = coco_add_event(opts, 'EP', 'boundary', par, '>', pint(2));
  end
end

%% construct finite state machine

states = { ...
  ... % initialisation states
  'init_prcond' 'init_chart' 'init_admissible' 'init_atlas' ...
  ... % covering algorithm
  'co_predict' 'co_correct' 'co_add_chart' 'co_flush' 'co_refine_step' ...
  ... % event handling
  'ev_init' ...
  'BP_locate' 'BP_add' 'BP_warning' ...
  'MX_check'  ...
  'SP_locate' 'SP_add' 'SP_warning' ...
  'ev_locate' 'ev_locate_cont' 'ev_locate_reg' ...
  'ev_locate_sing' 'ev_locate_multi' 'ev_warning' };

for i=1:numel(states)
	state                 = states{i};
	prop                  = state;
  full_state            = sprintf('state_%s', state);
  state_bslot           = sprintf('FSM_%s_begin', full_state);
  state_eslot           = sprintf('FSM_%s_end'  , full_state);
	han                   = str2func(full_state);
	cont.fsm.(prop).func  = han;
	cont.fsm.(prop).bslot = state_bslot;
	cont.fsm.(prop).eslot = state_eslot;
  opts = coco_add_signal(opts, state_bslot, 'covering');
  opts = coco_add_signal(opts, state_eslot, 'covering');
end
opts = coco_add_signal(opts, 'FSM_init', 'covering');

%% set entry point for continuer
cont.run  = @covering_run;

%% construct linear solver, corrector and atlas
[opts, cont, opts.atlas] = cont.atlas(opts, cont, dim);

% post merge defaults for corrector toolbox
cont = coco_merge(struct('corrector', 'nwtn'), cont);
if ischar(cont.corrector)
  cfname         = sprintf('corr_%s', cont.corrector);
  cont.corrector = str2func(cfname);
end
[opts, cont, opts.corr] = cont.corrector(opts, cont, cont.efunc);

% post merge defaults for linear solver toolbox
cont = coco_merge(struct('linsolve', 'splu'), cont);
if ischar(cont.linsolve)
  cfname        = sprintf('lsol_%s', cont.linsolve);
  cont.linsolve = str2func(cfname);
end
[opts, cont, opts.lsol] = cont.linsolve(opts, cont);

opts.cont = cont;

%% close system of equations and construct initial guess
[opts, opts.cont.chart0] = coco_close_efunc(opts, pars);
opts = coco_close_mfunc(opts);

%% initialise point list
opts = AtlasBase.bddat_init(opts);
opts = AtlasBase.bddat_set (opts, 'ins_mode', 'append');

end

%%
function [opts, cont] = get_settings(opts, dim)

defaults.NPR       = 10     ; % output every NPR steps
defaults.NSV       = []     ; % save solution every NSV steps, default = NPR
defaults.MEVFac    = 5      ; % tolerance factor for accepting multiple events

defaults.efunc     = 'efunc'; % function defining manifold
defaults.atlas     = []     ; % atlas toolbox
defaults.indcs     = opts.efunc.x_idx;

defaults.corr      = struct(); % post-init settings for corrector
defaults.lsol      = struct(); % post-init settings for linear solver

defaults.atlas_classes = { [] 'atlas_kd' ; 0 'atlas_0d' ; 1 'atlas_1d' };

cont = coco_get(opts, 'cont');
cont = coco_merge(defaults, cont);

assert(numel(cont.indcs) >= dim && max(cont.indcs) <= opts.efunc.x_dim, ...
  '%s: incompatible dimensions or variable indices', mfilename);

if isempty(cont.atlas)
  idx = find([cont.atlas_classes{:,1}]==dim,1) + 1;
  if isempty(idx)
    cont.atlas = str2func(sprintf('%s.create', cont.atlas_classes{  1,2}));
  else
    cont.atlas = str2func(sprintf('%s.create', cont.atlas_classes{idx,2}));
  end
elseif ischar(cont.atlas)
  cont.atlas = str2func(sprintf('atlas_%s.create', cont.atlas));
end

cont.filter = {'LogLevel' 'NPR' 'NSV' 'MEVFac' 'corr'};
end

%%
function opts = covering_run(opts, varargin)
%COCO_CONT  Entry point to continuation method.
%
%   OPTS = COCO_CONT(OPTS)
%   OPTS = COCO_CONT(OPTS,INITFLAG) The function COCO_CONT repeatedly calls
%   the function COCO_CONT_STEP which in turn executes the action of the
%   finite-state machine associated with the current state.  The COCO_CONT
%   function accepts (at least) two arguments, where the second argument
%   (when included) specifies whether this is a new run of the function or
%   an already initiated continuation. In the former case, the state of the
%   finite-state machine is set to 'init', whereas it is left unchanged in
%   the latter case.
%
%   To do: Add error handling for the case when compute_events detects
%   events that cannot be handled by handle_event or that should definitely
%   be addressed within a developer's toolbox. This will include adding an
%   'event' property to the 'cont' class.
%
%   See also: coco_cont_step
%

if nargin <= 1 || varargin{1}
	opts.cont.state = 'init_prcond';
end

opts.cont.accept = 0;

while(~opts.cont.accept)
	opts = covering_step(opts);
end

end

function [opts] = covering_step(opts)
%COCO_CONT_STEP  Dispatch call to current state in finite-state machine.

% execute state and emit 'begin' and 'end' signals
state = opts.cont.state;
fsm   = opts.cont.fsm;
opts  = coco_emit(opts, fsm.(state).bslot);
opts  = fsm.(state).func (opts);
opts  = coco_emit(opts, fsm.(state).eslot);

end
