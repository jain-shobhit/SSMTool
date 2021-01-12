function [opts chart f J] = efunc_FDF(opts, chart, xp)
%EFUNC_DFDX  Evaluate Jacobian of extended system at UP=[U;P].
%
%   [OPTS J] = EFUNC_DFDX(OPTS, UP) evaluates the Jacobian of the extended
%   system, which is the system combined of all algorithms and user
%   functions.
%
% See also: EFUNC_F
%

% create cache for variables that grow
persistent rows cols vals

if nargin==0
  % clear cache
  rows = [];
  cols = [];
  vals = [];
  return;
end

%% check if only evaluation of function required
if nargout<4
  [opts chart f] = efunc_F(opts, chart, xp);
  return
end

%% extract x and p from xp=[x;p] and
%  initialise f and
%  initialise derivatives with respect to x(efunc.p_idx,1)

efunc                = opts.efunc;
x                    = xp(efunc.xp_idx,1);
p                    = efunc.F_par_vals;
p(efunc.acp_f_idx,1) = x(efunc.p_idx,1);
f                    = efunc.f_init;
J                    = efunc.J_init;
[r c v]              = find(J);
[nrows ncols]        = size(J);
idx_MX               = numel(r);
idx                  = 1:idx_MX;
rows(idx)            = r;
cols(idx)            = c;
vals(idx)            = v;

%% compute Jacobian of extended continuation problem

for i = [ efunc.zero efunc.embedded ]
  
  func  = efunc.funcs(i);
  f_idx = func.f_idx;
  x_idx = func.x_idx;
  data  = func.data;
  
  % evaluate function and Jacobian of function
  switch func.call_mode
    case 1 % [d f]=F(o,d,x); [d J]=DFDX(o,d,x)
      [data F] = func.F(opts, data, x(x_idx));
      if isempty(func.DFDX)
        if func.vectorised
          [data DFDX] = coco_ezDFDX('f(o,d,x)v', opts, data, func.F, x(x_idx));
        else
          [data DFDX] = coco_ezDFDX('f(o,d,x)',  opts, data, func.F, x(x_idx));
        end
      else
        [data DFDX] = func.DFDX(opts, data, x(x_idx));
        if func.chkdrv
          if func.vectorised
            [data DFDX2] = coco_ezDFDX('f(o,d,x)v', opts, data, func.F, x(x_idx));
          else
            [data DFDX2] = coco_ezDFDX('f(o,d,x)',  opts, data, func.F, x(x_idx));
          end
          efunc_chkDeriv(opts, func, DFDX, DFDX2);
        end
      end
    case 2 % [d f J]=FDF(o,d,x)
      [data F DFDX] = func.F(opts, data, x(x_idx));
      if func.chkdrv
        if func.vectorised
          [data DFDX2] = coco_ezDFDX('f(o,d,x)v', opts, data, func.F, x(x_idx));
        else
          [data DFDX2] = coco_ezDFDX('f(o,d,x)',  opts, data, func.F, x(x_idx));
        end
        efunc_chkDeriv(opts, func, DFDX, DFDX2);
      end
    case 3 % [d c f]=F(o,d,c,x); [d c J]=DFDX(o,d,c,x)
      [data chart2 F] = func.F(opts, data, chart, x(x_idx));
      chart.private.data = chart2.private.data;
      if isempty(func.DFDX)
        if func.vectorised
          [data chart2 DFDX] = coco_ezDFDX('f(o,d,c,x)v', opts, data, chart, func.F, x(x_idx));
        else
          [data chart2 DFDX] = coco_ezDFDX('f(o,d,c,x)',  opts, data, chart, func.F, x(x_idx));
        end
      else
        [data chart2 DFDX] = func.DFDX(opts, data, chart, x(x_idx));
        if func.chkdrv
          if func.vectorised
            [data chart3 DFDX2] = coco_ezDFDX('f(o,d,c,x)v', opts, data, chart, func.F, x(x_idx)); %#ok<ASGLU>
          else
            [data chart3 DFDX2] = coco_ezDFDX('f(o,d,c,x)',  opts, data, chart, func.F, x(x_idx)); %#ok<ASGLU>
          end
          efunc_chkDeriv(opts, func, DFDX, DFDX2);
        end
      end
      chart.private.data = chart2.private.data;
    case 4 % [d c f J]=FDF(o,d,c,x)
      [data chart2 F DFDX] = func.F(opts, data, chart, x(x_idx));
      chart.private.data = chart2.private.data;
      if func.chkdrv
        if func.vectorised
          [data chart3 DFDX2] = coco_ezDFDX('f(o,d,c,x)v', opts, data, chart, func.F, x(x_idx)); %#ok<ASGLU>
        else
          [data chart3 DFDX2] = coco_ezDFDX('f(o,d,c,x)',  opts, data, chart, func.F, x(x_idx)); %#ok<ASGLU>
        end
        efunc_chkDeriv(opts, func, DFDX, DFDX2);
      end
  end
  
  if ~isempty(f_idx) % Matlab bug: this if should be unnecessary
    f(f_idx) = F - p(f_idx);
  end
  
  % merge Jacobians
  [r c v] = find(DFDX);
  idx       = idx_MX+(1:numel(r));
  idx_MX    = idx_MX+numel(r);
  rows(idx) = f_idx(r)';
  cols(idx) = x_idx(c)';
  vals(idx) = v;
  
  % allow modification of associated toolbox data
  % for example, counting evaluations of F and DFDX
  opts.efunc.funcs(i).data = data;
  
end

%% apply permutations
% [opts chart J1] = coco_ezDFDX('f(o,c,x)', opts, @efunc_F, chart, xp);

idx = 1:idx_MX;
f = f(efunc.f_idx);
J = sparse(efunc.f_idx(rows(idx)), efunc.xp_idx(cols(idx)), vals(idx), nrows, ncols);
