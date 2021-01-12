function [opts ev_vals] = efunc_events_F(opts, par_vals)
%EFUNC_EVENTS_F  Evaluate user functions (test functions).
%
%   [OPTS EVS] = EFUNC_EVENTS_F(OPTS, PARS)
%
% See also: EFUNC_F, EFUNC_DFDX, EFUNC_INIT
%

efunc = opts.efunc;

if ~isfield(efunc, 'ev') || isempty(efunc.ev.pidx)
  ev_vals = zeros(0,1);
	return;
end

ev          = efunc.ev;
ev_vals     = par_vals(ev.midx,1) - ev.vals;
