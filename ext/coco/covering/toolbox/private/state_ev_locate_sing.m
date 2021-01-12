function opts = state_ev_locate_sing(opts)
% use subdivision method to locate event on line segment

% initialise subdivision method

evidx = opts.cont.events.evidx;
cseg  = opts.cseg;

e0  = cseg.ev0;
e1  = cseg.ev1;

% [opts cseg chart] = cseg.chart_at(opts, 0, evidx);
% [opts e] = opts.efunc.events_F(opts, chart.p);
% if e(evidx)*e0<0
%   warning('possible spurious change');
% end

h0  = cseg.h0;
h1  = cseg.h1;

la1 = cseg.la1;
la2 = cseg.la2;

% compute first subdivision point

if abs(e0)<=abs(e1)
	h = la2*h0+la1*h1;
else
	h = la1*h0+la2*h1;
end

[opts cseg chart] = cseg.chart_at(opts, h, evidx);

%% locate special point

while norm(h1-h0)/(1+cseg.ptlist{1}.R) >= opts.cont.TOL
  [opts chart.e] = opts.efunc.events_F(opts, chart.p);
	e = chart.e(evidx);
	if e0*e<=0
		h1 = h;
		e1 = e;
	else
		h0 = h;
		e0 = e;
	end
	if e==0; break; end
	if abs(e0)<=abs(e1)
		h = la2*h0+la1*h1;
	else
		h = la1*h0+la2*h1;
  end
  [opts cseg chart] = cseg.chart_at(opts, h, evidx, chart);
end

% take right end point to guarantee that event occurs between
% cseg.u0 and chart.x (necessary for boundary events to avoid chopping off
% events at boundary)
[opts cseg chart] = cseg.chart_at(opts, h1, evidx, chart);
cseg.curr_chart   = chart;
opts.cseg         = cseg;

%% next state is add_<PointType>

opts.cont.state = opts.cont.locate_add_state;
