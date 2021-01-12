function opts = state_ev_init(opts)
% initialise curve segment for event handling

[opts opts.cseg]    = opts.cseg.ev_init(opts);
opts.cseg.src_chart = opts.cseg.ptlist{end};
opts.cont.events    = struct();
opts.cont.state     = 'BP_locate';
