function opts = state_init_prcond(opts)
% STATE_INIT_1  Initialise finite state machine, stage one.

% set global log priorities to initialise level
opts.cont.logPrioMN = opts.run.logPrioMN;
opts.cont.scrPrioMN = opts.run.scrPrioMN;
if opts.cont.LogLevel(1)>0
  opts.run.logPrioMN = opts.cont.LogLevel(1)+1;
  if numel(opts.cont.LogLevel)>=2
    opts.run.scrPrioMN = opts.cont.LogLevel(2)+1;
  else
    opts.run.scrPrioMN = opts.run.scrPrioMN+1;
  end
else
  opts.run.logPrioMN = 0;
  opts.run.scrPrioMN = 0;
end

% emit signal 'init'
opts = coco_emit(opts, 'FSM_init');

% initialise timer
opts.cont.tm = clock;

% reset corrector options to defaults
opts.cont.corr_opts = [];

% initialise first chart
[opts opts.atlas opts.cseg correct] = opts.atlas.init_prcond(opts, opts.cont.chart0);

% clear parametrised states
opts.cont.err_state  = [];
opts.cont.next_state = [];

% set next state
if correct
  % initialise corrector
  chart                 = opts.cseg.curr_chart;
  [opts chart accept x] = opts.corr.init(opts, chart, chart.x);
  chart.x               = x;
  opts.cseg.curr_chart  = chart;
  
  % set next state
  if accept
    opts.cont.state      = 'init_chart';
  else
    opts.cont.state      = 'co_correct';
    opts.cont.next_state = 'init_chart';
    opts.cont.err_state  = 'co_flush';
  end
else
  % accept point as is
  opts.cont.state = 'init_chart';
end
