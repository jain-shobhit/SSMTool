function [prob, chart, x0, V0, all_funcs] = coco_remesh(prob, chart, x0, V0, RMMX)

% return deep copy of functions if requested
if nargout>=5
  all_funcs = coco_save_funcs(prob); % not compatible with output syntax - Harry
end

if nargin<5
  RMMX = 10;
end

if isempty(V0)
  V0 = zeros(numel(x0),0);
end

% emit remesh signal to tell functions that they will be re-initialised very
% soon, for example, the projection condition shipped with CurveSegmentBase
prob = coco_emit(prob, 'remesh');

efunc = prob.efunc;
%% change all functions (and adjoints) and remesh x0 and V0

nefuncs  = numel(efunc.funcs);
old_x0  = x0(efunc.x_idx);
old_V0  = V0(efunc.x_idx,:);
old_x0p = x0(efunc.p_idx);
old_V0p = V0(efunc.p_idx,:);

efunc.chart = chart;
if ~isfield(prob, 'adjoint')
  for trial=1:RMMX
    efunc   = reset_efunc(efunc);
    status  = {};
    ids     = {};
    prob.efunc = efunc;
    for i=1:nefuncs
      prob.efunc.cfidx = i;
      func   = prob.efunc.funcs(i);
      [prob, s] = remesh_func(prob, func, x0, old_x0, old_V0, chart);
      if ~any(strcmp(s, {'success', 'fail', 'repeat'}))
        emsg = sprintf('%s: %s', mfilename, ...
          'remesh function returned illegal status:');
        error('%s\n%s returned with ''%s''', emsg, func.identifyer, s);
      end
      status = [ status, s ]; %#ok<AGROW>
      ids    = [ ids, func.identifyer ]; %#ok<AGROW>
    end
    efunc   = prob.efunc;
    success = strcmp('success', status);
    fail    = strcmp('fail'   , status);
    retry   = strcmp('repeat' , status);
    old_x0 = efunc.x0;
    old_V0 = efunc.V0;
    accept = all(success);
    if accept || any(fail)
      break
    end
  end
else
  adjoint = prob.adjoint;
  complementary = prob.complementary;
  old_l0  = old_x0(end-adjoint.a_dim(1)-complementary.v_dim+1:end);
  old_Vl0 = old_V0(end-adjoint.a_dim(1)-complementary.v_dim+1:end,:);
  for trial=1:RMMX
    efunc   = reset_efunc(efunc);
    adjoint = reset_adjoint(adjoint);
    complementary = reset_complementary(complementary);
    status  = {};
    ids     = {};
    prob.efunc   = efunc;
    prob.adjoint = adjoint;
    prob.complementary = complementary;
    for i=1:nefuncs
      prob.efunc.cfidx = i;
      func   = prob.efunc.funcs(i);
      [prob, s] = remesh_func(prob, func, x0, old_x0, old_V0, chart);
      if ~any(strcmp(s, {'success', 'fail', 'repeat'}))
        emsg = sprintf('%s: %s', mfilename, ...
          'remesh function returned illegal status:');
        error('%s\n%s returned with ''%s''', emsg, func.identifyer, s);
      end
      status = [ status, s ]; %#ok<AGROW>
      ids    = [ ids, func.identifyer ]; %#ok<AGROW>
      if i<=numel(adjoint.midx) && ~isempty(adjoint.midx{i})
        for j=adjoint.midx{i}
          prob.adjoint.cfidx = j;
          adjt = prob.adjoint.funcs(j);
          [prob, s] = remesh_adjt(prob, adjt, old_l0, old_Vl0, chart);
          if ~any(strcmp(s, {'success', 'fail', 'repeat'}))
            emsg = sprintf('%s: %s', mfilename, ...
              'adjoint remesh function returned illegal status:');
            error('%s\n%s returned with ''%s''', emsg, adjt.identifyer, s);
          end
          status = [ status, s ]; %#ok<AGROW>
          ids    = [ ids, adjt.identifyer ]; %#ok<AGROW>
        end
      end
    end
    efunc   = prob.efunc;
    adjoint = prob.adjoint;
    complementary = prob.complementary;
    success = strcmp('success', status);
    fail    = strcmp('fail'   , status);
    retry   = strcmp('repeat' , status);
    old_x0  = efunc.x0;
    old_V0  = efunc.V0;
    old_l0  = adjoint.l0;
    old_Vl0 = adjoint.Vl0;
    accept = all(success);
    if accept || any(fail)
      break
    end
  end
end

if ~accept
  emsg = sprintf('%s: remeshing failed:', mfilename);
  for i = find(fail | retry)
    emsg = sprintf('%s\n%s returned with status ''%s'' after %d trials', ...
      emsg, ids{i}, status{i}, RMMX);
  end
  error(emsg);
else
  coco_log(prob, 1, 1, '%s: mesh accepted after %d trials\n', ...
    mfilename, trial);
end
x0 = [ old_x0 ; old_x0p ];
V0 = [ old_V0 ; old_V0p ];

%% reclose efunc
cont_pars    = efunc.cont_pars;
cont_par_idx = coco_par2idx(prob, cont_pars, 'sloppy');

%% compute number of parameters to activate
% ppnum = m + d - n, m+d=efunc.f_dim, n = efunc.x_dim+efunc.p_dim
ppnum = efunc.f_dim-(efunc.x_dim+efunc.p_dim);

%% create copy of indices of potential continuation parameters
pararrays = { 'inactive_pars' 'active_pars' 'internal_pars' 'inequality_pars' };

for i = 1:numel(pararrays)
  tmp.(pararrays{i}) = efunc.(pararrays{i});
end

%% exchange parameters as defined with coco_xchg_pars
if isfield(efunc, 'xchg')
  xchg = coco_par2idx(prob, efunc.xchg);
  tmp  = exchange_pars(tmp, pararrays, xchg);
end

%% exchange internal parameters
pnum  = numel(cont_pars);         % number of parameters passed as argument
ipnum = numel(tmp.internal_pars); % number of internal pars
ipend = min(ipnum, pnum-ppnum);   % exchange internal parameters 1:ipend

if pnum>ppnum % more cont_pars than required -> exchange internal parameters
  xchg = [];
  for i = 1:ipend
    idx1 = tmp.internal_pars(i);
    idx2 = cont_par_idx(ppnum+i);
    xchg = [ xchg ; idx1 idx2 ]; %#ok<AGROW>
  end
  tmp = exchange_pars(tmp, pararrays, xchg);
end

%% compute subset of inactive parameters of argument cont_pars
inactive_cont_par_idx = [];
for idx = cont_par_idx
  if any(tmp.inactive_pars == idx)
    inactive_cont_par_idx = [ inactive_cont_par_idx idx ]; %#ok<AGROW>
  end
end
cpnum = numel(inactive_cont_par_idx); % number of inactive cont_pars

if cpnum<ppnum
  % we have more equations than variables+parameters
  errmsg = sprintf('%s: cannot close equations, too few parameters activated?\n', mfilename);
  errmsg = sprintf('%snumber of parameters to activate   : %d\n', errmsg, ppnum);
  errmsg = sprintf('%snumber of inactive parameters given: %d', errmsg, cpnum);
  error(errmsg); %#ok<SPERR>
end

%% activate ppnum inactive parameters

tmp.primary_pars = inactive_cont_par_idx(1:ppnum);
% for i = 1:ppnum
% 	% verify that cont_par_idx(1:ppnum) are now inactive
% 	if isempty(find(tmp.inactive_pars == cont_par_idx(i), 1))
% 		error('%s: too few inactive parameters specified for continuation', ...
% 			mfilename);
% 	end
% 	add parameter to list primary_pars
% 	tmp.primary_pars = [ tmp.primary_pars inactive_cont_par_idx(i)];
% end

efunc.acp_idx = [ ...
  inactive_cont_par_idx(1:ppnum) ...
  tmp.internal_pars ...
  tmp.active_pars ...
  tmp.inequality_pars];
efunc.acp_f_idx = efunc.pidx2fidx(efunc.acp_idx);

%% check for duplicate active continuation parameters

for i=1:numel(efunc.acp_idx)
  idx = find(efunc.acp_idx == efunc.acp_idx(i));
  if numel(idx) ~= 1
    error('%s: duplicate continuation parameter ''%s''', ...
      mfilename, efunc.idx2par{efunc.acp_idx(i)});
  end
end

%% create list of parameters for screen output

efunc.op_idx = [ cont_par_idx(1:ppnum) ...
  tmp.internal_pars(ipend+1:end) cont_par_idx(ppnum+1:end) ];

%% create permutation vector for xp

efunc.p_idx  = efunc.x_dim+(1:efunc.p_dim+ppnum);
efunc.xp_idx = [efunc.x_idx efunc.p_idx];
efunc.xp_dim = efunc.x_dim+efunc.p_dim+ppnum;

%% create event structure efunc.ev

ev.pidx       = [];
ev.midx       = [];
ev.par_type   = {};
ev.vals       = [];
ev.point_type = {};
ev.evgroup    = {};
ev.idx        = [];
ev.evsign     = '';
ev.BP_idx     = [];
ev.MX_idx     = [];
ev.SP_idx     = [];

% create copy of indices of remaining parameters
pararrays = { 'regular_pars' 'singular_pars' };
for i = 1:numel(pararrays)
  tmp.(pararrays{i}) = efunc.(pararrays{i});
end

% we ignore events in inactive parameters since these remain constant,
% but include primary continuation parameters
pararrays = { 'primary_pars' 'active_pars' 'inactive_pars' ...
  'internal_pars' 'inequality_pars' 'regular_pars' 'singular_pars' };
partypes  = { 'continuation' 'continuation' 'continuation'  ...
  'continuation' 'continuation' 'regular' 'singular' };

if isfield(efunc, 'events')
  events = efunc.events;
else
  events = [];
end

evnum = 0;

for i=1:numel(events)
  
  % look for parameters name in idx2par
  pars = events(i).par;
  pidx = [];
  for j=1:numel(pars)
    idx = find( strcmp(pars{j}, efunc.idx2par), 1 );
    if isempty(idx)
      error('%s: parameter ''%s'' not found, cannot add events', ...
        mfilename, pars{j});
    end
    pidx = [ pidx ; idx ]; %#ok<AGROW>
  end
  
  % compute type of event parameters
  ptype = {};
  for j=1:numel(pidx)
    type = {};
    for k=1:numel(pararrays)
      if ~isempty(find(tmp.(pararrays{k})==pidx(j),1))
        type = partypes{k};
        break
      end
    end
    if isempty(type)
      ptype = {};
      break;
    else
      ptype = [ ptype type ]; %#ok<AGROW>
    end
  end
  if isempty(ptype)
    continue
  end
  
  vals       = events(i).vals;
  vnum       = numel(vals);
  evidx      = evnum+(1:vnum);
  point_type = events(i).name;
  evgroup    = {[]};
  evsign     = events(i).sign;
  idx        = i;
  
  % expand arrays to match size of vals
  o            = ones(vnum, 1);
  point_type   = { point_type{o} };
  if numel(pidx)==1
    pidx       =   pidx(o);
    ptype      = { ptype{o}   };
    evsign     =   evsign(o)';
  else
    evgroup    = { evidx       };
  end
  evgroup      = { evgroup{o}  };
  idx          =   idx(o);
  
  % update entries in event structure
  ev.pidx       = [ ev.pidx         ; pidx                  ];
  ev.midx       = [ ev.midx         ; efunc.pidx2midx(pidx) ];
  ev.par_type   = [ ev.par_type       ptype                 ];
  ev.vals       = [ ev.vals         ; vals                  ];
  ev.point_type = [ ev.point_type     point_type            ];
  ev.evgroup    = [ ev.evgroup        evgroup               ];
  ev.evsign     = [ ev.evsign         evsign                ];
  ev.idx        = [ ev.idx          ; idx                   ];
  
  evlist      = events(i).evlist;
  ev.(evlist) = [ev.(evlist) evidx];
  
  evnum       = evnum + vnum;
end

efunc.ev   = ev;
prob.efunc = efunc;

% bug: remove argument chart and pass translation table ytr instead
[prob, chart] = efunc_remesh(prob, chart, x0);

end

%% local functions
function efunc = reset_efunc(efunc)

efunc.xtr       = [];
efunc.x0        = [];
efunc.V0        = [];
efunc.x_dim     =  0;
efunc.p_dim     =  0;
efunc.f_dim     =  0;
efunc.m_dim     =  0;
efunc.x_idx     = [];
efunc.f_idx     = [];
efunc.pidx2midx = [];
efunc.pidx2fidx = [];
efunc.idx2par   = {};
pararrays = { 'zero_pars' 'inactive_pars' 'active_pars' 'internal_pars' ...
  'inequality_pars' 'regular_pars' 'singular_pars' };
for i = 1:numel(pararrays)
  efunc.(pararrays{i}) = [];
end

end

function adjoint = reset_adjoint(adjoint)

adjoint.atr    = cell(1,2);
adjoint.a_dim  = zeros(1,2);
adjoint.ax_idx = [];
adjoint.af_idx = [];
adjoint.l0     = [];
adjoint.Vl0    = [];
adjoint.l_idx  = [];
adjoint.s_idx  = [];

end

function complementary = reset_complementary(complementary)

complementary.xtr     = [];
complementary.v0      = [];
complementary.V0      = [];
complementary.f_dim   =  0;
complementary.v_dim   =  0;
complementary.v_idx   = [];
complementary.f_idx   = [];
complementary.idx2par = {};

end

function [prob, s] = remesh_func(prob, func, x0, old_x0, old_V0, chart)

x0beg = numel(prob.efunc.x0);
xtr   = prob.efunc.xtr;

if isempty(func.remesh)
  if any(numel(func.x_idx) == numel(x0)) || strcmpi(func.x_idx, 'all')
    prob = coco_change_func(prob, func.data, 'xidx', 'all');
  else
    xidx = xtr(func.x_idx(func.x_idx<=numel(xtr)));
    if any(xidx<=0)
      error('%s: compulsory variables have been removed, cannot relocate ''%s''', ...
        mfilename, func.identifyer);
    end
    x0idx = func.x_idx(func.x_idx>numel(xtr));
    if isempty(xidx)
      prob = coco_change_func(prob, func.data, 'x0', old_x0(x0idx), ...
        'vecs', old_V0(x0idx,:));
    else
      prob = coco_change_func(prob, func.data, 'xidx', xidx, ...
        'x0', old_x0(x0idx), 'vecs', old_V0(x0idx,:));
    end
    xtr    = [ xtr ; x0beg+(1:numel(x0idx))' ];
  end
  s = 'success';
else
  x0idx = func.x_idx;
  [prob, s, tr] = func.remesh(prob, func.data, chart, old_x0(x0idx), ...
    old_V0(x0idx,:));
  xbeg  = x0beg*(tr~=0);
  xtr   = [ xtr ; xbeg(:)+tr(:) ];
end

prob.efunc.xtr = xtr;

end

function [prob, s] = remesh_adjt(prob, adjt, old_l0, old_Vl0, chart)

axbeg = numel(prob.adjoint.ax_idx);
afbeg = numel(prob.adjoint.af_idx);
[atr_x, atr_f] = deal(prob.adjoint.atr{:});
afidx = adjt.af_idx;
if isempty(adjt.remesh)
  aidx = atr_x(adjt.ax_idx(adjt.ax_idx<=numel(atr_x)));
  if any(aidx<=0)
    error('%s: compulsory variables have been removed, cannot relocate ''%s''', ...
      mfilename, adjt.identifyer);
  end
  axidx = adjt.ax_idx(adjt.ax_idx>numel(atr_x));
  if isempty(aidx)
    prob = coco_change_adjt(prob, adjt.data, 'l0', old_l0(afidx), ...
      'vecs', old_Vl0(afidx,:));
  else
    prob = coco_change_adjt(prob, adjt.data, 'aidx', aidx, ...
      'l0', old_l0(afidx), 'vecs', old_Vl0(afidx,:));
  end
  atr_x = [ atr_x ; axbeg+(1:numel(axidx))' ];
  atr_f = [ atr_f ; afbeg+(1:numel(afidx))' ];
  s = 'success';
else
  [prob, s, xtr, ftr] = adjt.remesh(prob, adjt.data, chart, old_l0(afidx), ...
    old_Vl0(afidx,:));
  xbeg  = axbeg*(xtr~=0);
  atr_x = [ atr_x ; xbeg(:)+xtr(:) ];
  fbeg  = afbeg*(ftr~=0);
  atr_f = [ atr_f ; fbeg(:)+ftr(:) ];
end

prob.adjoint.atr = { atr_x atr_f };

end

function tmp = exchange_pars(tmp, pararrays, xchg)

for i=1:size(xchg,1)
  
  for j=1:numel(pararrays)
    old_pararray = pararrays{j};
    old_pidx = find(tmp.(old_pararray)==xchg(i,1));
    if ~isempty(old_pidx); break; end
  end
  
  for j=1:numel(pararrays)
    new_pararray = pararrays{j};
    new_pidx = find(tmp.(new_pararray)==xchg(i,2));
    if ~isempty(new_pidx); break; end
  end
  
  pidx = tmp.(old_pararray)(old_pidx);
  tmp.(old_pararray)(old_pidx) = tmp.(new_pararray)(new_pidx);
  tmp.(new_pararray)(new_pidx) = pidx;
  
end

end
