function [prob cont corr] = corr_nwtns(prob, cont, func)

% set up damped Newton's method
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

%%
function [prob corr] = get_settings(prob)

defaults.ItMX       = 10    ; % max. number of iterations
defaults.SubItMX    = 4     ; % max. number of damping steps
defaults.ga0        = 1.0   ; % initial damping factor
defaults.al         = 0.5   ; % increase damping factor

corr                = coco_get(prob, 'corr');
corr                = coco_merge(defaults, corr);

corr.filter = {'LogLevel' 'TOL' 'ItMX' 'SubItMX' 'ga0' 'al' 'lsol'};

end

function corr = set_opts(corr, settings)
corr = coco_merge(corr, settings, corr.filter);
end

function settings = get_opts(corr)
settings = coco_merge(struct(), corr, corr.filter);
end

%%
function [prob chart x] = solve(prob, chart, x0)

[prob chart accept x] = init(prob, chart, x0);
while ~accept
	[prob chart accept x] = step(prob, chart);
end

end

%%
function [prob chart accept x0] = init(prob, chart, x0)

corr = prob.corr;

corr.It    = 0;
corr.SubIt = 0;
corr.ftm   = 0;
corr.dftm  = 0;
corr.stm   = 0;

prob.corr = corr;
accept    = false;

tm              = clock;
[prob chart f1] = corr.FDF(prob, chart, x0);
corr.ftm        = corr.ftm + etime(clock, tm);

prob = coco_emit(prob, 'corr_begin', 'nwtns', chart, x0);
prob = print_headline(prob, corr);
prob = print_data    (prob, corr, chart, x0, f1, 0, 0);

end

%%
function [prob chart accept x] = step(prob, chart)

corr = prob.corr;

corr.It = corr.It + 1;
x       = chart.x;

tm               = clock;
[prob chart f J] = corr.FDF(prob, chart, x);
corr.dftm        = corr.dftm + etime(clock, tm);

tm             = clock;
[prob chart d] = prob.lsol.solve(prob, chart, J, f);
corr.stm       = corr.stm + etime(clock, tm);

ga = corr.ga0;

for SubIt = 1:corr.SubItMX
	corr.SubIt      = SubIt;
	x               = chart.x - ga * d;
	tm              = clock;
	[prob chart f1] = corr.FDF(prob, chart, x);
	corr.ftm        = corr.ftm + etime(clock, tm);
  
	if (norm(f1)<=norm(f)) || (SubIt>=corr.SubItMX); break; end
  
	ga              = corr.al * ga;
end

prob.corr = corr;
accept    = (norm(d) < corr.TOL);

prob                 = print_data(prob, corr, chart, x, f1, d, ga);
[prob fids stop msg] = coco_emit(prob, 'corr_step', 'nwtns', chart, x);
stop                 = cell2mat(stop);

if accept
  prob = coco_emit(prob, 'corr_end', 'nwtns', 'accept', chart, x);
elseif ~isempty(stop) && any(stop)
  prob = coco_emit(prob, 'corr_end', 'nwtns', 'stop');
  emsg = sprintf('%s: stop requested by slot function(s)\n', mfilename);
  for idx=find(stop)
    emsg = sprintf('%s%s: %s\n', emsg, fids{idx}, msg{idx});
  end
	errmsg.identifier = 'CORR:Stop';
	errmsg.message = emsg;
	errmsg.FID = 'NWTNS';
	errmsg.ID  = 'MX';
	error(errmsg);
elseif corr.It >= corr.ItMX
  prob = coco_emit(prob, 'corr_end', 'nwtns', 'fail');
	errmsg.identifier = 'CORR:NoConvergence';
	errmsg.message = sprintf('%s %d %s', ...
		'no convergence of Newton''s method within', ...
		corr.ItMX, 'iterations');
	errmsg.FID = 'NWTNS';
	errmsg.ID  = 'MX';
	error(errmsg);
end

end

%%
function prob = print_headline(prob, corr) %#ok<INUSD>
%NWTN_PRINT_HEADLINE  Print headline for Newton's iteration.
%
%   OPTS = NWTN_PRINT_HEADLINE(OPTS) prints a headline for the iteration
%   information printed for each Newton step. This function calls
%   OPTS.NWTN.PRINT_HEADLINE to allow printing of a headline for
%   additional output data.
%
%   See also: print_data
%

coco_print(prob, 2, '\n%8s%10s%20s%10s%21s\n', ...
  'STEP', 'DAMPING', 'NORMS', ' ', 'COMPUTATION TIMES');
coco_print(prob, 2, '%4s%4s%10s%10s%10s%10s%7s%7s%7s', ...
  'IT', 'SIT', 'GAMMA', '||d||', '||f||', '||U||', 'F(x)', 'DF(x)', 'SOLVE');

prob = coco_emit(prob, 'corr_print', 'init', 2);
coco_print(prob, 2, '\n');

end

%%
function prob = print_data(prob, corr, chart, x, f, d, ga)
%NWTN_PRINT_DATA  Print information about Newton's iteration.
%
%   OPTS = NWTN_PRINT_DATA(OPTS) is called after each Newton iteration and
%   prints information its progress. This function calls
%   OPTS.NWTN.PRINT_DATA to allow printing of additional data.
%
%   See also: print_headline, coco_default_print
%

if corr.It==0
  coco_print(prob, 2, ...
    '%4d%4s%10s%10s%10.2e%10.2e%7.1f%7.1f%7.1f', ...
    corr.It, '', '', '', norm(f), norm(x), corr.ftm, corr.dftm, corr.stm);
else
  coco_print(prob, 2, ...
    '%4d%4d%10.2e%10.2e%10.2e%10.2e%7.1f%7.1f%7.1f', ...
    corr.It, corr.SubIt, ga, norm(d), norm(f), norm(x), ...
    corr.ftm, corr.dftm, corr.stm);
end

prob = coco_emit(prob, 'corr_print', 'data', 2, chart, x);
coco_print(prob, 2, '\n');

end
