function opts = state_co_predict(opts)
% STATE_PREDICT  Predict point and initialise corrector.

% create curve segment with predicted point in curr_chart
[opts opts.atlas opts.cseg correct] = opts.atlas.predict(opts, opts.cseg);

% set next state
if correct
  % initialise corrector
  chart                 = opts.cseg.curr_chart;
  [opts chart accept x] = opts.corr.init(opts, chart, chart.x);
  chart.x               = x;
  opts.cseg.curr_chart  = chart;
  
  if accept
    opts.cont.state      = 'co_add_chart';
    opts.cont.err_state  = [];
    opts.cont.next_state = [];
  else
    opts.cont.state      = 'co_correct';
    opts.cont.next_state = 'co_add_chart';
    opts.cont.err_state  = 'co_refine_step';
  end
else
  opts.cont.state        = 'co_add_chart';
end
