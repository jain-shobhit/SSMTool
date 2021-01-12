function opts = state_init_admissible(opts)
% Check if first chart is admissible. Flag admissible directions.

S             = struct('ep_flag', 0, 'ignore', [], 'msg', '', 'dir_flags', []);
cseg          = opts.cseg;
[opts cseg S] = check_point_and_directions(opts, cseg, S);
if S.ep_flag > 0
  cseg.curr_chart.ep_flag    = S.ep_flag;
  cseg.curr_chart.ignore_evs = S.ignore;
  [opts opts.atlas cseg]     = opts.atlas.init_admissible(opts, cseg, S);
end
opts.cseg = cseg;

if opts.cseg.Status == opts.cseg.CurveSegmentOK
  opts.cont.state = 'init_atlas';
else
  opts.cont.state = 'co_flush';
end

end

function [opts cseg S] = check_point_and_directions(opts, cseg, S)

evidx = union(opts.efunc.ev.BP_idx, opts.efunc.ev.MX_idx);
chart = cseg.curr_chart;

% evaluate event functions
cseg.curr_chart.t = nan(size(chart.x));
[opts cseg]       = cseg.eval_p(opts, evidx);
[opts e0]         = opts.efunc.events_F(opts, cseg.curr_chart.p);
if ~isfield(chart, 'p')
  chart.p = cseg.curr_chart.p;
end

% check point if at least one event was triggered
% we ignore monitor functions depending on t
% [this somewhat cumbersome argument list ensures that S.evidx is a row
% vector; this hack is necessary, because the set functions intersect,
% setdiff and union behave differently for different Matlab versions]
S.evidx = setdiff(evidx(:)', find(isnan(e0))');
S.e0    = e0;
if ~isempty(evidx)
  if isfield(chart, 'TS')
    % compute gradients of event functions
    h        = sqrt(opts.cont.TOL)*( 1.0 + max(abs(chart.x)) );
    x0       = chart.x;
    m        = size(chart.TS,2);
    grad_evF = zeros(numel(e0), m);
    for j = 1:m
      cseg.curr_chart.x = x0 + h*chart.TS(:,j);
      [opts cseg]       = cseg.eval_p(opts, evidx);
      [opts e1]         = opts.efunc.events_F(opts, cseg.curr_chart.p);
      grad_evF(:,j)     = (e1-e0)/h;
    end
    S.grad_evF        = grad_evF;
    cseg.curr_chart.t = chart.t;
    cseg.curr_chart.x = chart.x;
    [opts cseg S]     = check_point(opts, cseg, S);
    if S.ep_flag==1 && isfield(chart, 's') && ~isempty(chart.s)
      if isfield(chart, 'G')
        S.dfds_evF      = S.grad_evF * chart.G * chart.s;
      else
        S.dfds_evF      = S.grad_evF * chart.s;
      end
      [opts cseg S]   = check_directions(opts, cseg, S);
    end
  else
    cseg.curr_chart.t = chart.t;
    [opts cseg S]     = check_point(opts, cseg, S);
  end
else
  cseg.curr_chart.t = chart.t;
end

cseg.curr_chart.p = chart.p;

end

function [opts cseg S] = check_point(opts, cseg, S)

% compute subset of boundary and MX events that are triggered
evmask(opts.efunc.ev.BP_idx) = 1;
evmask(opts.efunc.ev.MX_idx) = 2;

evlist = intersect(find(evmask), S.evidx);
evsign = opts.efunc.ev.evsign(evlist);
evsidx = ((evsign=='<' | evsign=='>'))';
evlist = evlist(1,evsidx);
evsign = ((evsign=='<') - (evsign=='>'))';
evsign = evsign(evsidx,1);

epvals = S.e0(evlist);
N      = numel(evlist);
TOL    = ones(N,1)*opts.cont.TOL;
if isfield(S, 'grad_evF')
  dim = size(S.grad_evF,2);
  X   = eye(dim+1,dim);
  Y   = X;
  for i=1:N
    Y(dim+1,:) = S.grad_evF(evlist(i),:);
    TOL(i)     = (1+tan(subspace(X,Y)))*opts.cont.TOL;
  end
end
ephits = (evsign.*epvals <= TOL);
epidx  = find(ephits);
hitnum = numel(epidx);

% call event handler if defined to check event
for i=1:hitnum
  evidx = evlist(epidx(i));
  hanidx = opts.efunc.ev.idx(evidx);
  if ~isempty(opts.efunc.events(hanidx).han)
    msg.u      = chart.x;
    msg.e      = epvals(epidx(i));
    msg.pidx   = opts.efunc.ev.pidx(evidx);
    msg.midx   = opts.efunc.ev.midx(evidx);
    msg.p      = chart.p(msg.midx);
    msg.pars   = coco_idx2par(opts, msg.pidx);
    msg.evidx  = evidx;
    msg.action = '';
    
    evhan = opts.efunc.events(events.hanidx).han; % This looks like a typo: extra "events." here and on next three lines.
    data  = opts.efunc.events(events.hanidx).data;
    [data cseg events.msg] = evhan(opts, data, cseg, 'check_admissible', msg);
    opts.efunc.events(events.hanidx).data = data;
    
    if strcmp(msg.action, 'reject')
      ephits(epidx(i)) = false;
    end
  end
end

epidx    = find(ephits);
hitnum   = numel(epidx);
S.evidx  = evlist(epidx);
S.evsign = evsign(epidx);

% check if point is inside computational domain
if ~isempty(epidx)
  epnames = { opts.efunc.ev.point_type{evlist(epidx)} };
  msg = epnames{1};
  for i=2:hitnum
    msg = [ msg ', ' epnames{i} ]; %#ok<AGROW>
  end
  S.msg = msg;
  
  if all(evmask(evlist(ephits))==1) && all(abs(epvals(ephits))<=TOL(ephits))
    % chart on boundary
    S.ep_flag = 1;
    S.ignore  = evlist(ephits);
  else
    % chart outside boundary
    S.ep_flag = 2;
  end
end

end

function [opts cseg S] = check_directions(opts, cseg, S)

% compute subset of admissible directions
m      = size(S.dfds_evF,2);
flags  = ones(1,m)*opts.atlas.IsNotAdmissible;
evsign = S.evsign;

for i = 1:m
  ep = S.dfds_evF(S.evidx,i);
  if all(evsign.*ep >= 10*opts.cont.TOL)
    % direction points well inside computational domain within numerical
    % accuracy, a small step in this direction should always succeed
    flags(i) = opts.atlas.IsAdmissible;
  elseif all(evsign.*ep >= -10*opts.cont.TOL)
    % direction is within a possibly acceptable cone of near-tangent
    % directions, a step in this direction requires care
    flags(i) = opts.atlas.IsTangent;
  else
    % direction points well outside computational domain within numerical
    % accuracy, a step in this direction will usually fail
    flags(i) = opts.atlas.IsNotAdmissible;
  end
end

S.dir_flags = flags;

end
