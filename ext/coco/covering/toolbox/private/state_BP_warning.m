function opts = state_BP_warning(opts)

wmsg   = '';
events = opts.cont.events;
if isfield(events, 'hanidx')
  msg     = events.msg;
	pt_type = msg.point_type;
  if isfield(msg, 'wmsg')
    wmsg  = msg.wmsg;
  end
else
	pt_type = opts.efunc.ev.point_type{opts.cont.events.evidx};
end

if isempty(wmsg)
  coco_warn(opts, 1, opts.cont.LogLevel, ...
    '%s: error while locating boundary point\n', pt_type);
else
  coco_warn(opts, 1, opts.cont.LogLevel, ...
    '%s: error while locating boundary point', pt_type);
  coco_warn(opts, 1, opts.cont.LogLevel, '%s: %s\n', pt_type, wmsg);
end

if opts.cseg.ptlist{end}.ep_flag == 0
  opts.cseg.ptlist{end}.ep_flag = 2;
  opts.cseg.ptlist{end}.pt_type = 'MX';
end
opts.cseg.Status = opts.cseg.CurveSegmentOK;

opts.cont.events.state = 'init';
opts.cont.state        = 'ev_locate';
