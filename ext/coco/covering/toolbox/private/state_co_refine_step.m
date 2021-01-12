function opts = state_co_refine_step(opts)
%STATE ADD_POINT  Add regular point to curve segment.

[opts opts.atlas opts.cseg predict] = opts.atlas.refine_step(opts, opts.cseg);

if predict
  % reset MX flag to OK
  opts.cseg.Status = opts.cseg.CurveSegmentOK;
  
  opts.cont.state = 'co_predict';
else
  % next state is flush
  opts.cont.state = 'co_flush';
end
