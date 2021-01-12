function [prob cont corr] = corr_broyden(prob, cont, func)

% get toolbox settings
[prob corr] = get_settings(prob);

% set up Broyden's method
corr.F        = prob.(func).F;
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

corr = coco_get(prob, 'corr');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set defaults for Broyden's method

defaults.ItMX     = 25      ; % max. number of iterations
defaults.ItMN     = 3       ; % min. number of iterations
defaults.SubItMX  = 4       ; % max. number of damping steps
defaults.MaxStep  = 0.1     ; % max. size of Broyden step
defaults.ga0      = 1.0     ; % initial damping factor
defaults.al       = 0.5     ; % increase damping factor
defaults.lsol     = struct(); % settings for linear solver

corr = coco_merge(defaults, corr);
defaults.ResTOL = corr.TOL ; % convergence criterion ||f(x)||<=ResTOL
corr = coco_merge(defaults, corr);

corr.filter = {'LogLevel' 'TOL' 'ItMX' ...
  'ItMN' 'SubItMX' 'ResTOL' 'MaxStep' 'ga0' 'al' 'lsol'};

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
	[prob chart accept x] =step(prob, chart);
end

end

function [prob chart accept x0] = init(prob, chart, x0)

corr = prob.corr;

corr.x     = x0;
if ~isfield(corr, 'JInv')
  % bug: move JInv to chart data
  % bug: catch remesh event
  % bug: sparse updates (of factorisation?)
  corr.JInv  = eye(numel(x0));
end
corr.It    = 0;
corr.SubIt = 0;
corr.ftm   = 0;

tm                  = clock;
[prob chart corr.f] = corr.F(prob, chart, corr.x);
corr.ftm            = corr.ftm + etime(clock, tm);
corr.norm_f         = norm(corr.f);
corr.accept         = (corr.norm_f<prob.corr.ResTOL);
corr.accept         = corr.accept && (corr.It>=prob.corr.ItMN);

prob.corr = corr;
accept    = corr.accept;

prob = coco_emit(prob, 'corr_begin', 'broyden', chart, corr.x);
prob = print_headline(prob, corr);
prob = print_data    (prob, corr, chart, corr.x);
if accept
  prob = coco_emit(prob, 'corr_end', 'broyden', 'accept', chart, corr.x);
end

end

function [prob chart accept x] = step(prob, chart)

corr = prob.corr;

corr.It = corr.It + 1;

corr.d  = corr.JInv * corr.f;
corr.ga = corr.ga0;
x       = corr.x;
f       = corr.f;
if corr.ga*norm(corr.d) > corr.MaxStep*(1.0+norm(x))
	corr.ga = corr.MaxStep*(1.0+norm(x))/norm(corr.d);
end

for SubIt = 1:corr.SubItMX
	corr.SubIt          = SubIt;
	corr.x              = x - corr.ga * corr.d;
	tm                  = clock;
	[prob chart corr.f] = corr.F(prob, chart, corr.x);
	corr.ftm            = corr.ftm + etime(clock, tm);
	if (norm(corr.f)<max(corr.ResTOL, corr.norm_f)) || (SubIt>=corr.SubItMX)
		break;
	end
	corr.ga = corr.al * corr.ga;
end

corr.norm_f = norm(corr.f);
corr.accept = (norm(corr.norm_f) < corr.ResTOL);
corr.accept = corr.accept || (norm(corr.d) < corr.TOL);
corr.accept = corr.accept && (corr.It>=prob.corr.ItMN);

if ~corr.accept
  s = corr.x - x;
  y = corr.f - f;
  B = corr.JInv;
  r = (s'*B*y);
  if abs(r)>=corr.TOL
    corr.JInv = B + (s-B*y) * (s'*B) / r;
  end
end

prob.corr = corr;
x         = corr.x;
accept    = corr.accept;

prob                 = print_data(prob, corr, chart, corr.x);
[prob fids stop msg] = coco_emit(prob, 'corr_step', 'broyden', chart, corr.x);
stop                 = cell2mat(stop);

if accept
  prob = coco_emit(prob, 'corr_end', 'broyden', 'accept', chart, corr.x);
elseif ~isempty(stop) && any(stop)
  prob = coco_emit(prob, 'corr_end', 'broyden', 'stop');
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
  prob = coco_emit(prob, 'corr_end', 'broyden', 'fail');
	errmsg.identifier = 'CORR:NoConvergence';
	errmsg.message = sprintf('%s %d %s', ...
		'no convergence of Broyden''s method within', ...
		corr.ItMX, 'iterations');
	errmsg.FID = 'NWTN';
	errmsg.ID  = 'MX';
	error(errmsg);
end

end

function prob = print_headline(prob, corr) %#ok<INUSD>
%NWTN_PRINT_HEADLINE  Print headline for Newton's iteration.
%
%   OPTS = NWTN_PRINT_HEADLINE(OPTS) prints a headline for the iteration
%   information printed for each Newton step. This function calls
%   OPTS.NWTN.PRINT_HEADLINE to allow printing of a headline for
%   additional output data.
%
%   See also: broyden_print_data
%

coco_print(prob, 2, '\n%8s%10s%20s%10s%7s\n', ...
  'STEP', 'DAMPING', 'NORMS', ' ', 'TIME');
coco_print(prob, 2, ...
  '%4s%4s%10s%10s%10s%10s%7s', 'IT', 'SIT', ...
	'GAMMA', '||d||', '||f||', '||U||', 'F(x)');

prob = coco_emit(prob, 'corr_print', 'init', 2);
coco_print(prob, 2, '\n');

end

function prob = print_data(prob, corr, chart, x)
%NWTN_PRINT_DATA  Print information about Newton's iteration.
%
%   OPTS = NWTN_PRINT_DATA(OPTS) is called after each Newton iteration and
%   prints information its progress. This function calls
%   OPTS.NWTN.PRINT_DATA to allow printing of additional data.
%
%   See also: broyden_print_headline, coco_default_print
%

if corr.It==0
  coco_print(prob, 2, '%4d%4s%10s%10s%10.2e%10.2e%7.1f', ...
    corr.It, '', '', '', norm(corr.f), norm(corr.x), corr.ftm);
else
  coco_print(prob, 2, '%4d%4d%10.2e%10.2e%10.2e%10.2e%7.1f', ...
    corr.It, corr.SubIt, corr.ga, norm(corr.d), norm(corr.f), norm(corr.x), corr.ftm);
end

prob = coco_emit(prob, 'corr_print', 'data', 2, chart, x);
coco_print(prob, 2, '\n');

end
