function opts = state_ev_locate_multi(opts)
% use subdivision method to locate event on line segment

cseg = opts.cseg;

% loop over all indices

x = [];

evnum = numel(opts.cont.events.evidx);
for i=1:evnum
  
  % initialise subdivision method
  evidx = opts.cont.events.evidx(i);
  
  e0  = cseg.ev0(i);
  e1  = cseg.ev1(i);
  
  % ignore "touching zeros" and constant event functions
  if e0*e1>0 || abs(e1-e0)<=10*opts.cont.TOL
    continue
  end
  
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
  
  % locate special point
  
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
  [opts cseg chart] = cseg.chart_at(opts, h1, opts.cont.events.evidx, chart);
  [opts chart.e]    = opts.efunc.events_F(opts, chart.p);
  
  % end loop over all indices
  x(:,end+1) = chart.x; %#ok<AGROW>
end

cseg.curr_chart = chart;
opts.cseg = cseg;

% check if points are close to each other

xnum = size(x,2);

if xnum==0
  opts.cont.state = opts.cont.locate_warn_state;
  return;
  
elseif xnum==1
  xx  = x;
  err = 0;
  
else
  xx  = sum(x,2)/xnum;
  dx  = x-xx(:,ones(1,xnum));
  err = sqrt(max(sum(dx.*dx,1)));
  
end

% check if solution is acceptable and go to next state

if err <= 10*opts.cont.TOL
  evs = norm(chart.e(opts.cont.events.evidx));
  if evs <= 10*opts.cont.TOL
    chart.x         = xx;
    cseg.curr_chart = chart;
    opts.cseg       = cseg;
    opts.cont.state = opts.cont.locate_add_state;
  else
    opts.cont.state = opts.cont.locate_warn_state;
  end
else
  opts.cont.state = opts.cont.locate_warn_state;
end
