function opts = state_init_atlas(opts)
% STATE_INIT_3  Initialise finite state machine, stage 3.

% initialise atlas code with first complete chart
[opts opts.atlas opts.cseg flush] = opts.atlas.init_atlas(opts, opts.cseg);

% set next state
if opts.cseg.Status == opts.cseg.CurveSegmentOK
  if flush
    opts.cont.state = 'co_flush';
  else
    opts.cseg       = [];
    opts.cont.state = 'co_predict';
  end
else
  opts.cont.state = 'co_flush';
end
