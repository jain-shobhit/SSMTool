function opts = state_init_chart(opts)
% STATE_INIT_2  Initialise finite state machine, stage 2.

% update tangent space of first chart and initialise chart.s
[opts opts.atlas opts.cseg] = opts.atlas.init_chart(opts, opts.cseg);

% set post-initialise logging priorities
if opts.cont.LogLevel(1)>0
  opts.run.logPrioMN = opts.cont.LogLevel(1);
  if numel(opts.cont.LogLevel)>=2
    opts.run.scrPrioMN = opts.cont.LogLevel(2);
  else
    opts.run.scrPrioMN = opts.cont.scrPrioMN;
  end
end

% set post-initialise toolbox settings
opts.corr = opts.corr.set_opts(opts.corr, opts.cont.corr);
opts.lsol = opts.lsol.set_opts(opts.lsol, opts.cont.lsol);

% clear parametrised states
opts.cont.err_state  = [];
opts.cont.next_state = [];

% set next state
if opts.cseg.Status == opts.cseg.CurveSegmentOK
  if isa(opts.cseg, 'EventLocator')
    opts.cont.state = 'init_admissible';
  else
    opts.cont.state = 'init_atlas';
  end
else
  opts.cont.state = 'co_flush';
end
