function [prob cont corr] = corr_continex(prob, cont, func)

% set up Newton's method
[prob corr] = get_settings(prob);

corr.FDF      = prob.(func).FDF;
corr.solve    = @solve;
corr.init     = @init;
corr.step     = @step;
corr.set_opts = @set_opts;
corr.get_opts = @get_opts;

prob = coco_add_signal(prob, 'corr_begin', mfilename);
prob = coco_add_signal(prob, 'corr_step', mfilename);
prob = coco_add_signal(prob, 'corr_end', mfilename);

prob = coco_add_signal(prob, 'corr_print', mfilename);

prob = coco_add_signal(prob, 'corr_sample', mfilename);
prob = coco_add_signal(prob, 'corr_pull_back', mfilename, '-max-slots', 1);

end

function [prob corr] = get_settings(prob)

defaults.ItMX       = 10      ; % max. number of iterations
defaults.ItMN       = 0       ; % min. number of iterations
defaults.ItNW       = []      ; % max. number of full Newton iterations
defaults.SubItMX    = 4       ; % max. number of damping steps
defaults.NPt        = 3       ; % number of accepted points to avarage over
defaults.MaxStep    = 0.1     ; % max. relative size of Newton step
defaults.MaxAbsStep = inf     ; % max. absolute size of Newton step
defaults.MaxAbsDist = inf     ; % max. diameter of trust region
defaults.NIt        = 1       ; % number of iterations for computing abs step
defaults.MADSteps   = 0       ; % minimum number of steps to cross trust region
defaults.DampRes    = []      ; % use damping if ||f(x_It)||>DampRes(It)
defaults.ga0        = 1.0     ; % initial damping factor
defaults.al         = 0.5     ; % increase damping factor
defaults.sfac       = 0.95    ; % safety scaling factor
defaults.phan       = []      ; % axes handles for debug plots

corr                = coco_get(prob, 'corr');
corr                = coco_merge(defaults, corr);

defaults.ResTOL     = corr.TOL;   % convergence criteria:
defaults.corrMX     = corr.TOL;   %   (||f(x)||<=ResTOL && ||d||<=corrMX)
defaults.ResMX      = corr.TOL;   %   (||d||<=TOL && ||f(x)||<=ResMX)
defaults.varMX      = corr.TOL/2; % radius around mean

corr                = coco_merge(defaults, corr);

defaults.DampResMX  = corr.ResTOL; % use damping if ||f(x)||>DampResMX

corr                = coco_merge(defaults, corr);

corr.NIt      = max(1, corr.NIt);
corr.MADSteps = min(corr.MADSteps, corr.ItMX);

corr.filter = {'LogLevel' 'TOL' 'ItMX' ...
  'ItMN' 'ItNW' 'SubItMX' 'NPt' 'ResTOL' 'corrMX' 'varMX' 'ResMX' 'MaxStep' ...
  'MaxAbsStep' 'MaxAbsDist' 'NIt' 'MADSteps' 'DampRes' 'DampResMX' 'ga0' 'al' ...
  'sfac' 'lsol'};

end

function corr = set_opts(corr, settings)
corr = coco_merge(corr, settings, corr.filter);
end

function settings = get_opts(corr)
settings = coco_merge(struct(), corr, corr.filter);
end

function [prob chart x] = solve(prob, chart, x0)

[prob chart accept x] = init(prob, chart, x0);
while ~accept
	[prob chart accept x] = step(prob, chart);
end

end

function [prob chart accept x0] = init(prob, chart, x0)

corr = prob.corr;

corr.x      = x0;
corr.pts    = { x0 };
corr.It     = 0;
corr.ItPT   = ceil(linspace(corr.ItMN+(corr.ItMX-corr.ItMN)/corr.NPt,corr.ItMX,corr.NPt));
corr.SubIt  = 0;
corr.ftm    = 0;
corr.dftm   = 0;
corr.stm    = 0;
corr.sample = false;

tm                  = clock;
[prob chart corr.f] = corr.FDF(prob, chart, corr.x);
corr.ftm            = corr.ftm + etime(clock, tm);
corr.norm_f_old     = norm(corr.f);
corr.accept         = (corr.norm_f_old<prob.corr.ResTOL);
corr.accept         = corr.accept && (corr.It>=corr.ItMN);

if corr.accept
  corr.Pt        = 1;
  corr.X         = x0;
  corr.F         = corr.f;
  corr.accept    = (corr.Pt>=corr.NPt);
  corr.print.mrk = sprintf('%d', corr.Pt);
  corr.sample    = true;
  prob           = coco_emit(prob, 'corr_sample', corr.x);
else
  corr.Pt        = 0;
  corr.X         = [];
  corr.F         = [];
  corr.print.mrk = ' ';
end

prob.corr = corr;
accept    = corr.accept;

prob = coco_emit(prob, 'corr_begin', 'nwtn', chart, x0);
prob = print_headline(prob, corr);
prob = print_data    (prob, corr, chart, x0, 0, 0);
if accept
  prob = coco_emit(prob, 'corr_end', 'nwtn', 'accept', chart, x0);
end

coco_log(prob, 1, corr.LogLevel, ...
  '%s: MaxAbsStep = %.2e, MaxAbsDist = %.2e\n', ...
  mfilename, corr.MaxAbsStep, corr.MaxAbsDist);

end

function [prob chart accept x] = step(prob, chart)

corr = prob.corr;

corr.It = corr.It + 1;

if isempty(corr.ItNW) || corr.It <= corr.ItNW
	tm                         = clock;
	[prob chart corr.f corr.J] = corr.FDF(prob, chart, corr.x);
	corr.dftm                  = corr.dftm + etime(clock, tm);
else
	tm                  = clock;
	[prob chart corr.f] = corr.FDF(prob, chart, corr.x);
	corr.ftm            = corr.ftm + etime(clock, tm);
end

tm                  = clock;
[prob chart corr.d] = prob.lsol.solve(prob, chart, corr.J, corr.f);
corr.stm            = corr.stm + etime(clock, tm);

debug_plot(prob, corr, corr.f, corr.x, corr.d);

corr.ga  = corr.ga0;
x        = corr.x;
scale    = 1.0 ./ (1.0+abs(x));
rel_step = norm(scale.*corr.d)/(1+norm(scale.*x));
if corr.ga*rel_step > corr.MaxStep
	corr.ga = corr.MaxStep/rel_step;
end

[m n] = size(corr.X);
beg   = max(1, numel(corr.pts)+1-corr.NIt);
pts   = [ corr.pts(beg:end) mat2cell(corr.X, m, repmat(1,1,n)) ];
xx    = x - corr.ga*corr.d;
[prob fids u] = coco_emit(prob, 'corr_pull_back', xx); %#ok<ASGLU>
if numel(u)==1
  xx = u{1};
end
func = @(z) norm(z - xx);
dist = min(cellfun(func, pts));
while dist > corr.MaxAbsStep
	corr.ga = corr.sfac*corr.ga*corr.MaxAbsStep/dist;
  xx      = x - corr.ga*corr.d;
  [prob fids u] = coco_emit(prob, 'corr_pull_back', xx); %#ok<ASGLU>
  if numel(u)==1
    xx = u{1};
  end
  func = @(z) norm(z - xx);
  dist = min(cellfun(func, pts));
end

xx = x - corr.ga*corr.d;
[prob fids u] = coco_emit(prob, 'corr_pull_back', xx); %#ok<ASGLU>
if numel(u)==1
  xx = u{1};
end
dist = norm(corr.pts{1} - xx);
SubIt = 1;
while dist > corr.MaxAbsDist
	SubIt=SubIt+1;
	corr.ga = corr.al * corr.ga;
	if SubIt>corr.SubItMX
    prob = coco_emit(prob, 'corr_end', 'nwtn', 'fail');
    errmsg.identifier = 'CORR:NoConvergence';
    errmsg.message = sprintf('correction leaves trust reagion');
    errmsg.FID = 'NWTN';
    errmsg.ID  = 'MX';
    error(errmsg);
	end
  xx = x - corr.ga*corr.d;
  [prob fids u] = coco_emit(prob, 'corr_pull_back', xx); %#ok<ASGLU>
  if numel(u)==1
    xx = u{1};
  end
  dist = norm(corr.pts{1} - xx);
end

if numel(corr.DampRes)<corr.It
  DampResMX = corr.DampResMX;
else
  DampResMX = corr.DampRes(corr.It);
end

for SubIt = SubIt:corr.SubItMX
  corr.SubIt          = SubIt;
  dx                  = corr.ga * corr.d;
  corr.x              = x - dx;
  [prob fids u]       = coco_emit(prob, 'corr_pull_back', corr.x); %#ok<ASGLU>
  if numel(u)==1
    dx                = u{1}-x;
    corr.x            = u{1};
  end
  tm                  = clock;
  [prob chart corr.f] = corr.FDF(prob, chart, corr.x);
  corr.ftm            = corr.ftm + etime(clock, tm);
  norm_av_f           = norm(mean([corr.F corr.f], 2));
  accept              = ...
    ( (norm(corr.d) < corr.TOL) && (norm_av_f <= corr.ResMX) );
  accept              = accept || ...
    ( (norm_av_f < corr.ResTOL) && (norm(corr.d) <= corr.corrMX) );
  accept              = accept || ...
    (norm_av_f <= max(corr.norm_f_old, DampResMX));
  if accept || (SubIt>=corr.SubItMX)
    break;
  end
  corr.ga = corr.al * corr.ga;
end

x               = corr.x;
corr.pts        = [corr.pts x];
corr.norm_f_old = norm_av_f;
corr.accept     = ...
  ( (norm(corr.d) < corr.TOL) && (corr.norm_f_old <= corr.ResMX) );
corr.accept     = corr.accept || ...
  ( (corr.norm_f_old < corr.ResTOL) && (norm(corr.d) <= corr.corrMX) );
corr.accept     = corr.accept && ...
  (corr.It>=corr.ItMN);

if corr.accept
  if ~corr.sample
    corr.sample = true;
    prob        = coco_emit(prob, 'corr_sample', corr.x);
  end
  corr.Pt = corr.Pt+1;
  corr.X  = [ corr.X x ];
  corr.F  = [ corr.F corr.f ];
  dist    = corr.X - repmat(mean(corr.X, 2), 1, size(corr.X,2));
  dist    = sqrt(sum(dist.^2,1));
  Pts     = sum(dist<=corr.varMX);
  corr.accept = (Pts>=corr.NPt);
  if corr.accept
    x      = mean(corr.X, 2);
    corr.x = x;
  end
  corr.print.mrk = sprintf('%d', Pts);
else
  corr.print.mrk = ' ';
end

corr_MX = corr.It+max(corr.NPt-corr.Pt,1-corr.accept) > corr.ItMX;
corr_MX = corr_MX || corr.It >= corr.ItPT(min(corr.Pt+1,corr.NPt));
if corr_MX
  corr.print.mrk = 'MX';
  corr.accept    = false;
end

prob.corr = corr;
accept    = corr.accept;

prob                 = print_data(prob, corr, chart, x, norm(dx), norm(corr.pts{1} - x));
[prob fids stop msg] = coco_emit(prob, 'corr_step', 'nwtn', chart, x);
stop                 = cell2mat(stop);

if accept
  prob = coco_emit(prob, 'corr_end', 'nwtn', 'accept', chart, x);
elseif ~isempty(stop) && any(stop)
  prob = coco_emit(prob, 'corr_end', 'nwtn', 'stop');
  emsg = sprintf('%s: stop requested by slot function(s)\n', mfilename);
  for idx=find(stop)
    emsg = sprintf('%s%s: %s\n', emsg, fids{idx}, msg{idx});
  end
	errmsg.identifier = 'CORR:Stop';
	errmsg.message = emsg;
	errmsg.FID = 'NWTN';
	errmsg.ID  = 'MX';
	error(errmsg);
elseif corr_MX
  prob = coco_emit(prob, 'corr_end', 'nwtn', 'fail');
	errmsg.identifier = 'CORR:NoConvergence';
	errmsg.message = sprintf('%s %d %s', ...
		'no convergence of Newton''s method within', ...
		corr.ItMX, 'iterations');
	errmsg.FID = 'NWTN';
	errmsg.ID  = 'MX';
	error(errmsg);
end

end

function prob = print_headline(prob, corr) %#ok<INUSD>
coco_print(prob, 2, '\n%8s%3s%10s%20s%10s%21s\n', ...
  'STEP', ' ', 'DAMPING', 'NORMS', ' ', 'COMPUTATION TIMES');
coco_print(prob, 2, '%4s%4s%3s%10s%10s%10s%10s%7s%7s%7s%10s%10s', ...
  'IT', 'SIT', '#', 'GAMMA', '||d||', '||f||', '||U||', 'F(x)', 'DF(x)', 'SOLVE', '||dx||', 'dist');

prob = coco_emit(prob, 'corr_print', 'init', 2);
coco_print(prob, 2, '\n');
end

function prob = print_data(prob, corr, chart, x, dx, dst)
if corr.It==0
  coco_print(prob, 2, ...
    '%4d%4s%3s%10s%10s%10.2e%10.2e%7.1f%7.1f%7.1f%10s%10s', ...
    corr.It, '', corr.print.mrk, '', '', corr.norm_f_old, norm(corr.x), ...
    corr.ftm, corr.dftm, corr.stm, '', '');
else
  if isnan(dx)
    coco_print(prob, 2, ...
      '%4d%4d%3s%10.2e%10.2e%10.2e%10.2e%7.1f%7.1f%7.1f%10s%10s', ...
      corr.It, corr.SubIt, corr.print.mrk, corr.ga, ...
      norm(corr.d), corr.norm_f_old, norm(corr.x), ...
      corr.ftm, corr.dftm, corr.stm, '', '');
  else
    coco_print(prob, 2, ...
      '%4d%4d%3s%10.2e%10.2e%10.2e%10.2e%7.1f%7.1f%7.1f%10.2e%10.2e', ...
      corr.It, corr.SubIt, corr.print.mrk, corr.ga, ...
      norm(corr.d), corr.norm_f_old, norm(corr.x), ...
      corr.ftm, corr.dftm, corr.stm, dx, dst);
  end
end

prob = coco_emit(prob, 'corr_print', 'data', 2, chart, x);
coco_print(prob, 2, '\n');

end

function debug_plot(prob, corr, f, x, d)

switch numel(corr.phan)
  case 0
  case 1
    coco_plot_F(corr.phan(1), prob, f, 'b.-');
  otherwise
    coco_plot_F(corr.phan(1), prob, f, 'b.-');
    coco_plot_u(corr.phan(2), prob, x, 'b.-');
    hold(corr.phan(2), 'on')
    coco_plot_u(corr.phan(2), prob, -d, 'r.-');
    hold(corr.phan(2), 'off')
end
drawnow

end
