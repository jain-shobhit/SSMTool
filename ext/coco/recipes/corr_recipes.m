function [prob cont corr] = corr_recipes(prob, cont, func)
%CORR_RECIPES   Recipes for continuation: nonlinear solver.
%
% Set up the Newton-based nonlinear solver by assigning content to
% the corr return argument. Let the associated linear solver be
% given by lsol_recipes.m by default. Declare the 'corr_begin',
% 'corr_stop', 'corr_end', and 'corr_print' signals.
%
% Required subfunctions:
% get_settings - Return solver settings in corr argument.
% set_opts     - Merge selected settings.
% get_opts     - Extract selected settings.
% init         - Trivial encoding.
% solve        - Solve a nonlinear system using an iterated Newton scheme.
% step         - Take a single Newton step.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: corr_recipes.m 2839 2015-03-05 17:09:01Z fschild $

[prob corr] = get_settings(prob);
defaults.linsolve  = 'recipes' ; % linear solver toolbox
cont = coco_merge(defaults, cont);

corr.FDF      = prob.(func).FDF;
corr.solve    = @solve;
corr.init     = @init;
corr.step     = @step;
corr.set_opts = @set_opts;
corr.get_opts = @get_opts;

prob = coco_add_signal(prob, 'corr_begin', mfilename);
prob = coco_add_signal(prob, 'corr_step', mfilename);
prob = coco_add_signal(prob, 'corr_end', mfilename);

prob = coco_add_signal(prob, 'corr_print', 'corr_nwtns (recipes)');

end

%% Settings
function [prob corr] = get_settings(prob)

defaults.ItMX       = 10    ; % maximum number of iterations
defaults.SubItMX    = 8     ; % maximum number of damping steps
defaults.ga0        = 1.0   ; % initial damping factor
defaults.al         = 0.5   ; % increase in damping factor

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

%% Solve
function [prob chart x] = solve(prob, chart, x0)
%SOLVE   Iterated calls to step until accept = true.

[prob chart accept x] = init(prob, chart, x0);
while ~accept
	[prob chart accept x] = step(prob, chart);
end

end

function [prob chart accept x0] = init(prob, chart, x0)
%INIT   Initialize solver.
%
% The function initializes counters and timers, emits the
% 'corr_begin' signal, and prints initial data to screen.

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

function [prob chart accept x] = step(prob, chart)
%STEP   Take a single damped Newton step.
%
% Rescale the damping factor at most SubItMX times. Let acccept =
% true if the norm of the correction is below the global tolerance
% TOL. Emits the 'corr_step' signal. If a stop message is received
% or the number of iterates exceeds ItMX, then emits the 'corr_end'
% signal with associated error message.

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
  
	if norm(f1) <= norm(f); break; end
  
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

%% Screen output
function prob = print_headline(prob, corr) %#ok<INUSD>
%PRINT_HEADLINE  Print headline for Newton's iteration.
%
% The function emits the 'corr_print' signal in order to allow
% printing of a headline for additional output data.
%
% See also: print_data

coco_print(prob, 2, '\n%8s%10s%20s%10s%21s\n', ...
  'STEP', 'DAMPING', 'NORMS', ' ', 'COMPUTATION TIMES');
coco_print(prob, 2, '%4s%4s%10s%10s%10s%10s%7s%7s%7s', ...
  'IT', 'SIT', 'GAMMA', '||d||', '||f||', '||U||', 'F(x)', 'DF(x)', 'SOLVE');

prob = coco_emit(prob, 'corr_print', 'init', 2);
coco_print(prob, 2, '\n');

end

function prob = print_data(prob, corr, chart, x, f, d, ga)
%PRINT_DATA   Print information about progress of each Newton iteration.
%
% The function emits the 'corr_print' signal in order to allow
% printing of additional data.
%
% See also: print_headline, coco_print

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
