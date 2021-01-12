function opts = efunc_minit(opts)
%EFUNC_INIT  Initialise extended system.
%
%   OPTS = EFUNC_INIT(OPTS, X)
%
% See also: EFUNC_F, EFUNC_DFDX, EFUNC_MONITOR_F
%

%% affected fields in opts
%
%    opts.efunc.F_par_vals - values of user-parameters
%    opts.efunc.f_init     - 
%    opts.efunc.m_init     - 
%    opts.efunc.J_init     - 

efunc = opts.efunc;

%% initialise preallocated vectors of function values
%  for extended continuation problem and monitor function

efunc.m_init = nan(efunc.m_dim, 1);

%% postinitialise x-indices if required

xp_idx = 1:efunc.xp_dim;
for i = 1:numel(efunc.funcs)
  if strcmpi(efunc.funcs(i).x_idx, 'all')
    efunc.funcs(i).x_idx = xp_idx;
  end
end

%% copy efunc back to opts

opts.efunc = efunc;

