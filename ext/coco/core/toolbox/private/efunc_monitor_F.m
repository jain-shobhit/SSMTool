function [opts chart par_vals] = efunc_monitor_F(opts, chart, xp, txp, ev_idx)
%EFUNC_USER_F  Evaluate monitor functions (test functions).
%
%   [OPTS PARS] = EFUNC_MONITOR_F(OPTS, X)
%
% See also: EFUNC_F, EFUNC_DFDX, EFUNC_INIT
%

%% compute current values of monitor functions

efunc    = opts.efunc;
x        = xp (efunc.xp_idx);
t        = txp(efunc.xp_idx);
par_vals = efunc.m_init;
funcs    = [ efunc.embedded efunc.regular efunc.singular ];

if nargin<5
  for i = funcs
    func  = efunc.funcs(i);
    m_idx = func.m_idx;
    x_idx = func.x_idx;
    data  = func.data;
    
    % evaluate function
    [opts data chart F]  = efunc_call_F(opts, data, chart, func, x(x_idx), t(x_idx));
    par_vals(m_idx) = F;
    
    % allow modification of associated toolbox data
    % for example, counting evaluations of F
    opts.efunc.funcs(i).data = data;
  end
else
  ev_midx = efunc.ev.midx(ev_idx);
  for i = funcs
    func  = efunc.funcs(i);
    m_idx = func.m_idx;
    x_idx = func.x_idx;
    data  = func.data;
    
    if ~isempty(intersect(func.req_idx, ev_midx))
      % evaluate function
      [opts data chart F]  = efunc_call_F(opts, data, chart, func, x(x_idx), t(x_idx));
      par_vals(m_idx) = F;
      
      % allow modification of associated toolbox data
      % for example, counting evaluations of F
      opts.efunc.funcs(i).data = data;
    end
  end
end
