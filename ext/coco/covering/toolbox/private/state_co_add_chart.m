function opts = state_co_add_chart(opts)
%STATE ADD_POINT  Add regular point to atlas segment.

% add chart to atlas
[opts opts.atlas opts.cseg flush] = opts.atlas.add_chart(opts, opts.cseg);

% go to next state
if flush
  if opts.cseg.Status==opts.cseg.CurveSegmentOK ...
      && isa(opts.cseg, 'EventLocator')
    opts.cont.state = 'ev_init';
  else
    opts.cont.state  = 'co_flush';
  end
else
  opts.cont.state = 'co_predict';
end
