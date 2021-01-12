function [opts, chart] = efunc_init(opts)
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

efunc.f_init = nan(efunc.f_dim, 1);

%% postinitialise x-indices if required

xp_idx = 1:efunc.xp_dim;
for i = 1:numel(efunc.funcs)
  if strcmpi(efunc.funcs(i).x_idx, 'all')
    efunc.funcs(i).x_idx = xp_idx;
  end
end

%% initialise xp0, tp0 and nontrivial parameter values

% get explicitly defined parameter values (see coco_set_parival)
if isfield(efunc, 'parivals') && ~isempty(efunc.parivals)
  ipars =   efunc.parivals(:,1)   ;
  ivals = [ efunc.parivals{:,2} ]';
  ipidx = efunc.pidx2fidx(coco_par2idx(opts, ipars));
else
  ivals = [];
  ipidx = [];
end

F_par_vals = zeros(efunc.f_dim, 1);
x0         = efunc.x0;
chart = efunc.chart;

% compute initial values of embedded monitor functions
opts.efunc = efunc;
for i = efunc.embedded
  func = efunc.funcs(i);
  
  % compute initial parameter values
  f_idx             = func.f_idx;
  x_idx             = func.x_idx;
  data              = func.data;
  [opts, data, chart, pval] = efunc_call_F(opts, data, chart, func, x0(x_idx), []);
  F_par_vals(f_idx) = pval;
  
  % allow modification of associated toolbox data
  % for example, counting evaluations of F
  efunc.funcs(i).data = data;
  opts.efunc          = efunc;
  
end
F_par_vals(ipidx) = ivals;
xp0  = [ efunc.x0 ; F_par_vals(efunc.acp_f_idx) ];
tp0  = [ efunc.tx ; efunc.tp(efunc.acp_f_idx)   ];
if isfield(efunc, 'TSX')
  TS0 = [ efunc.TSX ; efunc.TSP(efunc.acp_f_idx,:) ];
end
efunc.F_par_vals = F_par_vals;

%% apply permutations and construct initial chart

xp0([efunc.x_idx efunc.p_idx],1) = xp0;
tp0([efunc.x_idx efunc.p_idx],1) = tp0;
if isfield(efunc, 'TSX')
  TS0([efunc.x_idx efunc.p_idx],:) = TS0;
end

chart.x = xp0;
chart.t = tp0;
if isfield(efunc, 'TSX')
  chart.TS = TS0;
end

%% construct numerical differentiation functions
%  bug: this needs to be filled in

%% initialise Jacobian of efunc with parameter derivatives

J_init       = speye(efunc.f_dim, efunc.f_dim);
efunc.J_init = [sparse(efunc.f_dim, efunc.x_dim) -J_init(:,efunc.acp_f_idx)];

%% copy efunc back to opts

opts.efunc = efunc;
