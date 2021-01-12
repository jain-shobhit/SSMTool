function opts = state_ev_locate_cont(opts)
% locate event using Newton's method on extended system

% initialise xfunc (see also compute_evlist)

ev     = opts.efunc.ev;
cseg   = opts.cseg;

evidx  = opts.cont.events.evidx;
mf_idx = find(opts.efunc.acp_idx==ev.pidx(evidx),1);
mf_idx = opts.efunc.p_idx(mf_idx);
mf_val = ev.vals(evidx);

% initialise nwtn

e0 = cseg.ev0;
e1 = cseg.ev1;
% u0 = events.u0;
% u1 = events.u1;
[opts cseg chart hh]  = cseg.chart_at(opts, -e0/(e1-e0), evidx);
opts                  = coco_emit(opts, 'update_h', hh);
opts                  = coco_emit(opts, 'fix_mfunc', mf_idx, mf_val);
[opts chart accept x] = opts.corr.init(opts, chart, chart.x);
chart.x               = x;
cseg.curr_chart       = chart;
opts.cseg             = cseg;

if accept
  % add point
  opts.cont.state      = opts.cont.locate_add_state;
else
  % correct, then add point
  opts.cont.state      = 'co_correct';
  opts.cont.next_state = opts.cont.locate_add_state;
  opts.cont.err_state  = opts.cont.locate_warn_state;
end
