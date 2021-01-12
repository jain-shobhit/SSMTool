function opts = state_co_flush(opts)
%STATE_CO_FLUSH  Flush point list.

% set post-initialise logging priorities, must be repeated here, because
% flush might be entered after init_prcond if correct fails
if opts.cont.LogLevel(1)>0
  opts.run.logPrioMN = opts.cont.LogLevel(1);
  if numel(opts.cont.LogLevel)>=2
    opts.run.scrPrioMN = opts.cont.LogLevel(2);
  else
    opts.run.scrPrioMN = opts.cont.scrPrioMN;
  end
end

% ignore future events that occur exactly at points in ptlist to avoid
% duplicate detections
cseg = opts.cseg;
m    = numel(cseg.ptlist);
if cseg.Status == cseg.CurveSegmentOK && ~cseg.isInitialSegment
  for i=1:m
    cseg.ptlist{i}.ignore_at = true;
  end
end
opts.cseg = cseg;

% flush point list of current curve segment
[opts opts.atlas opts.cseg] = opts.atlas.flush(opts, opts.cseg);

% stop FSM if cseg.Status not equal CurveSegmentOK
opts.cont.accept = (opts.cseg.Status~=opts.cseg.CurveSegmentOK);

% cleanup: restore log priorities before exiting FSM
if opts.cont.accept
  opts.run.logPrioMN = opts.cont.logPrioMN;
  opts.run.scrPrioMN = opts.cont.scrPrioMN;
end

% instance of cseg becomes invalid after flush -> destroy
% opts.atlas.base_chart = opts.cseg.base_chart;
opts.cseg = [];

% next state is predict
opts.cont.state = 'co_predict';

end
