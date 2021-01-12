function opts = state_co_correct(opts)
% STATE_CORRECT  Execute correction steps.

%% perform one newton iteration
try
  chart                 = opts.cseg.curr_chart;
  [opts chart accept x] = opts.corr.step(opts, chart);
  chart.x               = x;
  opts.cseg.curr_chart  = chart;
catch err
  % mark curve segment as invalid and handle error
  opts.cseg.Status = opts.cseg.CorrectionFailed;
	switch err.identifier
    case 'CORR:NoConvergence'
      if isempty(opts.cont.err_state)
        rethrow(err);
      else
        coco_warn(opts, 3, opts.cont.LogLevel, ...
          '%s: correction failed\n%s\n', mfilename, err.message);
        opts.cont.state      = opts.cont.err_state;
        opts.cont.err_state  = [];
        opts.cont.next_state = [];
        return
      end
      
		otherwise
      rethrow(err);
	end
end

if accept
	opts.cont.state      = opts.cont.next_state;
	opts.cont.err_state  = [];
	opts.cont.next_state = [];
end
