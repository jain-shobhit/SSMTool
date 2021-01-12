function opts = state_ev_locate_reg(opts)
% locate event using hybrid subdivision + Newton method

% evaluate event function at new subdivision point

evidx          = opts.cont.events.evidx;
cseg           = opts.cseg;
[opts cseg]    = cseg.eval_p(opts, evidx);
chart          = cseg.curr_chart;
[opts chart.e] = opts.efunc.events_F(opts, chart.p);

v   = chart.x;
h   = cseg.h;
e0  = cseg.ev0;
evs = chart.e;
e   = evs(evidx);

% compute new subdivision interval

if e0*e<=0
	cseg.v1  = v;
	cseg.h1  = h;
	cseg.ev1 = e;
else
	cseg.v0  = v;
	cseg.h0  = h;
	cseg.ev0 = e;
end

% check convergence

if abs(cseg.h1-cseg.h0)/(1+cseg.ptlist{1}.R) <= opts.cont.TOL || e==0
  % take right end point to guarantee that event occurs between
  % cseg.u0 and chart.x (necessary for boundary events to avoid chopping off
  % events at boundary)
  chart.x         = cseg.v1;
  % bug: comment the next 4 lines
  [opts chart]         = opts.cseg.update_TS(opts, chart);
  [opts chart]         = opts.cseg.update_t(opts, chart);
  [opts chart chart.p] = opts.efunc.monitor_F(opts, chart, chart.x, chart.t);
  [opts chart.e]       = opts.efunc.events_F (opts, chart.p);
  % [opts cseg chart] = cseg.chart_at(opts, cseg.h1, evidx, chart); % Harry added

  cseg.curr_chart = chart;
  opts.cseg       = cseg;
	% go to add_<PointType>
	opts.cont.state = opts.cont.locate_add_state;
	return
end

% initialise subdivision method
e0 = cseg.ev0;
e1 = cseg.ev1;

v0 = cseg.v0;
h0 = cseg.h0;
v1 = cseg.v1;
h1 = cseg.h1;

la1 = cseg.la1;
la2 = cseg.la2;

% compute subdivision point
if abs(e0)<=abs(e1)
	v = la2*v0+la1*v1;
	h = la2*h0+la1*h1;
else
	v = la1*v0+la2*v1;
	h = la1*h0+la2*h1;
end
cseg.h = h;

% interpolate
[opts cseg chart hh] = cseg.chart_at(opts, h, evidx, chart);

% initialise xfunc
opts = coco_emit(opts, 'update_h', hh);

% initialise nwtn
[opts chart f1] = opts.efunc.F(opts, chart, chart.x);
[opts chart f2] = opts.efunc.F(opts, chart, v);

if norm(f1)<norm(f2)
  [opts chart accept x] = opts.corr.init(opts, chart, chart.x);
else
  [opts chart accept x] = opts.corr.init(opts, chart, v);
end

chart.x         = x;
cseg.curr_chart = chart;
opts.cseg       = cseg;

% go to next state
if accept
  % go to locate_reg
  opts.cont.state      = 'ev_locate_reg';
else
  % next state is correct, then go to locate_reg
  opts.cont.state      = 'co_correct';
  opts.cont.next_state = 'ev_locate_reg';
  opts.cont.err_state  = opts.cont.locate_warn_state;
end
