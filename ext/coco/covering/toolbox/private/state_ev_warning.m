function opts = state_ev_warning(opts)

coco_warn(opts, 2, opts.cont.LogLevel, ...
  'error while locating event of type ''%s''\n', ...
  opts.efunc.ev.point_type{evidx});

opts.cont.state = 'ev_locate';
