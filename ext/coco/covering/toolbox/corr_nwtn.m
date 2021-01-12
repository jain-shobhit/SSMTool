function [prob cont corr] = corr_nwtn(prob, cont, func)

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

end

function [prob corr] = get_settings(prob)

defaults.ItMX       = 10      ; % max. number of iterations
defaults.ItMN       = 0       ; % min. number of iterations
defaults.ItNW       = []      ; % max. number of full Newton iterations
defaults.SubItMX    = 4       ; % max. number of damping steps
defaults.MaxStep    = 0.1     ; % max. relative size of Newton step
defaults.MaxAbsStep = inf     ; % max. absolute size of Newton step
defaults.MaxAbsDist = inf     ; % max. diameter of trust region
defaults.DampRes    = []      ; % use damping if ||f(x_It)||>DampRes(It)
defaults.ga0        = 1.0     ; % initial damping factor
defaults.al         = 0.5     ; % increase damping factor
defaults.phan       = []      ; % axes handles for debug plots

corr                = coco_get(prob, 'corr');
corr                = coco_merge(defaults, corr);

defaults.ResTOL     = corr.TOL; % convergence criteria:
defaults.corrMX     = corr.TOL; %   (||f(x)||<=ResTOL && ||d||<=corrMX)
defaults.ResMX      = corr.TOL; %   (||d||<=TOL && ||f(x)||<=ResMX)

corr                = coco_merge(defaults, corr);

defaults.DampResMX  = corr.ResTOL; % use damping if ||f(x)||>DampResMX

corr                = coco_merge(defaults, corr);

corr.filter = {'LogLevel' 'TOL' 'ItMX' ...
  'ItMN' 'ItNW' 'SubItMX' 'ResTOL' 'corrMX' 'ResMX' 'MaxStep' ...
  'MaxAbsStep' 'MaxAbsDist' 'DampRes' 'DampResMX' 'ga0' 'al' 'lsol'};

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

corr.x     = x0;
corr.pts   = { x0 };
corr.It    = 0;
corr.SubIt = 0;
corr.ftm   = 0;
corr.dftm  = 0;
corr.stm   = 0;

tm                  = clock;
[prob chart corr.f] = corr.FDF(prob, chart, corr.x);
corr.ftm            = corr.ftm + etime(clock, tm);
corr.norm_f_old     = norm(corr.f);
corr.accept         = (corr.norm_f_old<prob.corr.ResTOL);
corr.accept         = corr.accept && (corr.It>=corr.ItMN);

prob.corr = corr;
accept    = corr.accept;

prob = coco_emit(prob, 'corr_begin', 'nwtn', chart, x0);
prob = print_headline(prob, corr);
prob = print_data    (prob, corr, chart, x0);
if accept
  prob = coco_emit(prob, 'corr_end', 'nwtn', 'accept', chart, x0);
end

end

function [prob chart accept x] = step(prob, chart)

  function dist = func1(pts,x,ga,d)
    out = zeros(numel(pts),1);
    for i=1:numel(pts)
      out(i) = norm(pts{i} - x + ga*d);
    end
    dist = min(out);
  end

  function out = func2(ga,pt,x,d)
    out = norm(pt - x + ga*d);
  end
    
corr = prob.corr;

corr.It = corr.It + 1;

if isempty(corr.ItNW) || corr.It <= corr.ItNW
	tm                         = clock;
	[prob chart corr.f corr.J] = corr.FDF(prob, chart, corr.x);
	corr.dftm                  = corr.dftm + etime(clock, tm);
  % J2 = coco_ezDFDX('f(x,p)', @(x,p) debug_F(prob, corr, chart, x), corr.x, nan(0,1));
else
	tm                  = clock;
	[prob chart corr.f] = corr.FDF(prob, chart, corr.x);
	corr.ftm            = corr.ftm + etime(clock, tm);
end

tm                  = clock;
[prob chart corr.d] = prob.lsol.solve(prob, chart, corr.J, corr.f);
corr.stm            = corr.stm + etime(clock, tm);

if ~isempty(corr.phan)
  debug_plot(prob, corr, corr.f, corr.x, corr.d);
end

corr.ga  = corr.ga0;
x        = corr.x;
scale    = 1.0 ./ (1.0+abs(x));
rel_step = norm(scale.*corr.d)/(1+norm(scale.*x));
if corr.ga*rel_step > corr.MaxStep
	corr.ga = corr.MaxStep/rel_step;
end

% func = @(z) norm(z - x + corr.ga*corr.d);
% dist = min(cellfun(func, corr.pts));
if func1(corr.pts, x, corr.ga, corr.d) > corr.MaxAbsStep
	corr.ga = corr.MaxAbsStep/(norm(corr.d));
end
% 
%func  = @(ga) norm(corr.pts{1} - x + ga*corr.d);
SubIt = 1;
while func2(corr.ga, corr.pts{1}, x, corr.d) > corr.MaxAbsDist
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
end

if numel(corr.DampRes)<corr.It
  DampResMX = corr.DampResMX;
else
  DampResMX = corr.DampRes(corr.It);
end

for SubIt = SubIt:corr.SubItMX
	corr.SubIt          = SubIt;
	corr.x              = x - corr.ga * corr.d;
	tm                  = clock;
	[prob chart corr.f] = corr.FDF(prob, chart, corr.x);
	corr.ftm            = corr.ftm + etime(clock, tm);
  accept              = ...
    ( (norm(corr.d) < corr.TOL) && (norm(corr.f) <= corr.ResMX) );
  accept              = accept || ...
    ( (norm(corr.f) < corr.ResTOL) && (norm(corr.d) <= corr.corrMX) );
  accept              = accept || ...
      (norm(corr.f) <= max(corr.norm_f_old, DampResMX));
	if accept || (SubIt>=corr.SubItMX)
		break;
  end
	corr.ga = corr.al * corr.ga;
end

x               = corr.x;
corr.pts        = [corr.pts x];
corr.norm_f_old = norm(corr.f);
corr.accept     = ...
  ( (norm(corr.d) < corr.TOL) && (corr.norm_f_old <= corr.ResMX) );
corr.accept     = corr.accept || ...
  ( (corr.norm_f_old < corr.ResTOL) && (norm(corr.d) <= corr.corrMX) );
corr.accept     = corr.accept && ...
  (corr.It>=corr.ItMN);

prob.corr = corr;
accept    = corr.accept;

prob                 = print_data(prob, corr, chart, x);
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
elseif corr.It >= corr.ItMX
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
coco_print(prob, 2, '\n%8s%10s%20s%10s%21s\n', ...
  'STEP', 'DAMPING', 'NORMS', ' ', 'COMPUTATION TIMES');
coco_print(prob, 2, '%4s%4s%10s%10s%10s%10s%7s%7s%7s', ...
  'IT', 'SIT', 'GAMMA', '||d||', '||f||', '||U||', 'F(x)', 'DF(x)', 'SOLVE');

prob = coco_emit(prob, 'corr_print', 'init', 2);
coco_print(prob, 2, '\n');
end

function prob = print_data(prob, corr, chart, x)
if corr.It==0
  coco_print(prob, 2, ...
    '%4d%4s%10s%10s%10.2e%10.2e%7.1f%7.1f%7.1f', ...
    corr.It, '', '', '', norm(corr.f), norm(corr.x), ...
    corr.ftm, corr.dftm, corr.stm);
else
  coco_print(prob, 2, ...
    '%4d%4d%10.2e%10.2e%10.2e%10.2e%7.1f%7.1f%7.1f', ...
    corr.It, corr.SubIt, corr.ga, ...
    norm(corr.d), norm(corr.f), norm(corr.x), ...
    corr.ftm, corr.dftm, corr.stm);
end

prob = coco_emit(prob, 'corr_print', 'data', 2, chart, x);
coco_print(prob, 2, '\n');
end

function f = debug_F(prob, corr, chart, x) %#ok<DEFNU>
	[prob chart f] = corr.FDF(prob, chart, x); %#ok<ASGLU>
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
