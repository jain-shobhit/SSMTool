function [opts chart f] = efunc_F(opts, chart, xp)
%EFUNC_F  Evaluate extended system at UP=[U;P].
%
%   [OPTS F] = EFUNC_F(OPTS, UP) evaluates the extended system, which is
%   the system combined of all algorithms and user functions.
%
% See also: EFUNC_DFDX
%

%% extract x and p from xp=[x;p] and initialise f

efunc                = opts.efunc;
x                    = xp(efunc.xp_idx,1);
p                    = efunc.F_par_vals;
p(efunc.acp_f_idx,1) = x(efunc.p_idx,1);
f                    = efunc.f_init;

%% evaluate extended continuation problem

for i = [ efunc.zero efunc.embedded ]
  
  func  = efunc.funcs(i);
  f_idx = func.f_idx;
  x_idx = func.x_idx;
  data  = func.data;
  
  % evaluate function
  [opts data chart F] = efunc_call_F(opts, data, chart, func, x(x_idx), []);
  
  if ~isempty(f_idx) % Matlab bug: this if should be unnecessary
    f(f_idx) = F - p(f_idx);
  end
  
  % allow modification of associated toolbox data
  % for example, counting evaluations of F
  opts.efunc.funcs(i).data = data;
  
end

%% apply permutation

f = f(efunc.f_idx);
