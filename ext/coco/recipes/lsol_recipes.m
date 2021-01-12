function [prob cont lsol] = lsol_recipes(prob, cont)
%LSOL_RECIPES   Recipes for continuation: linear solver.
%
% Set up the linear solver by assigning content to the lsol return
% argument. If either of the det or cond fields is true, then add
% chart data with identifier 'lsol' to each chart. If the cond
% field is true, append a slot function that prints the condition
% number to screen.
%
% Required subfunctions:
% get_settings - Return solver settings in lsol argument.
% set_opts     - Merge selected settings.
% get_opts     - Extract selected settings.
% init         - Trivial encoding.
% solve        - Solve a linear system using Matlab's backslash operator.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: lsol_recipes.m 2839 2015-03-05 17:09:01Z fschild $

[prob lsol] = get_settings(prob);

lsol.init     = @init;
lsol.solve    = @solve;
lsol.set_opts = @set_opts;
lsol.get_opts = @get_opts;

if lsol.det || lsol.cond
  prob = coco_add_chart_data(prob, 'lsol', [], []);
  if lsol.cond
    prob = coco_add_slot(prob, 'lsol.cond', @corr_print, [], 'corr_print');
  end
  if lsol.det
    sd = coco_func_data(struct('save_det', false));
    prob = coco_add_slot(prob, 'lsol.det', @set_save_det, sd, 'set_mode');
    lsol.sd = sd;
  end
end

end

%% Settings
function [prob lsol] = get_settings(prob)

defaults.thresh   = [0.1 0.001]; % pivot threshholds for UMFPACK
defaults.NFac     = 1;           % number of factors for scaled determinant
defaults.warn     = true;        % display warning
defaults.det      = false;       % monitor determinant in chart data
defaults.cond     = false;       % monitor condition number in chart data

lsol = coco_get(prob, 'lsol');
lsol = coco_merge(defaults, lsol);

lsol.filter = {'LogLevel' 'thresh' 'NFac' 'warn'};

end

function lsol = set_opts(lsol, settings)
lsol = coco_merge(lsol, settings, lsol.filter);
end

function settings = get_opts(lsol)
settings = coco_merge(struct(), lsol, lsol.filter);
end

%% Solver
function [prob chart] = init(prob, chart)
end

function [prob chart x] = solve(prob, chart, A, b)
%SOLVE  Solve A*x=b using Matlab's backslash operator.
%
% Linear solver (see 'doc mldivide' for further explanation).
% Optionally appends rescaled determinant and condition number of A
% to the det and cond chart data fields associated with identifier
% 'lsol'. For det close to zero, the v chart data field contains a
% unit vector perpendicular to the corresponding tangent vector.

S = warning('query', 'backtrace');
warning('off', 'backtrace');
[L,U,P,Q,R] = lu(A, prob.lsol.thresh);
if isempty(b)
  x = [];
else
  x = Q*(U\(L\(P*(R\b))));
end
warning(S.state, 'backtrace');

% check matrix format
[m n]    = size(A);
issquare = (m==n);

if prob.lsol.warn && ~issquare
  coco_warn(prob, 3, prob.lsol.LogLevel, ...
    'lsol: matrix A not square, size(A) = [%d %d]\n%s\n', m, n, ...
    'computation of determinant and condition number disabled');
  prob.lsol.warn = false;
end

cdata = struct();

% compute rescaled determinant
if prob.lsol.det && prob.lsol.sd.save_det
  if issquare
    DR  = full(diag(R));
    SR  = prod(sign(DR));
    DU  = full(diag(U));
    SU  = prod(sign(DU));
    DPQ = det(P)*det(Q);
    sru = sort(sort(abs(DR)).*sort(abs(DU),'descend'));
    
    % det(A) = DPQ*SR*SU*prod(sru)
    cdata.fdet = DPQ*SR*SU*prod(sru);
    % we rescale and multiply only a few of the smallest factors of prod(sru)
    N   = min(prob.lsol.NFac, m);
    la  = sru(1:N);
    la  = sqrt(la.*flipud(la));
    sc = abs(la);
    la = (2*la)./(max(1,sc)+sc);
    % sc  = nthroot(N,N/mean(sc))*(sc./( 1+sc ));
    % sc  = max(sru(1:N), ones(N,1))+sru(1:N);
    % sc  = 1+sru(1:N);
    % sc  = 2*( sru(1:N)./sc );
    cdata.det = DPQ*SR*SU*prod(la);
    [v1 idx1] = min(abs(DR));
    [v2 idx2] = min(abs(DU));
    if v1<=v2
      zb       = zeros(size(P,1),1);
      zb(idx1) = 1;
      v        = Q*(U\(L\(P*zb)));
    else
      UU            = U;
      UU(idx2,idx2) = 1;
      zb            = zeros(size(P,1),1);
      zb(idx2)      = 1;
      v             = Q*(UU\zb);
    end
    cdata.v = v/norm(v);
    % fprintf(2, 'det(J) = % .2e, scdet(A) = % .2e\n', ...
    %   DPQ*SR*SU*prod(sru), DPQ*SR*SU*prod(sc));
  else
    cdata.det = nan;
  end
end

if prob.lsol.cond
  if issquare
    ad = abs(diag(R));
    cn = max(ad)/min(ad);
    ad = abs(diag(L));
    cn = cn * max(ad)/min(ad);
    ad = abs(diag(U));
    cn = cn * max(ad)/min(ad);
    cn = full(cn);
  else
    cn = nan;
  end
  if prob.lsol.cond
    cdata.cond = cn;
  end
end

if prob.lsol.det || prob.lsol.cond
  chart = coco_set_chart_data(chart, 'lsol', cdata);
end

end

%% Screen output
function prob = corr_print(prob, data, cmd, LogLevel, chart, varargin) %#ok<INUSL>
%CORR_PRINT   Slot function that responds to the corr_print signal.

switch cmd
  case 'init'
    coco_print(prob, LogLevel, '%10s', 'cond(DF)');
  case 'data'
    cdata = coco_get_chart_data(chart, 'lsol');
    if isfield(cdata, 'cond')
      coco_print(prob, LogLevel, '%10.2e', cdata.cond);
    else
      coco_print(prob, LogLevel, '%10s', '');
    end
end
end

function data = set_save_det(prob, data, mode, varargin)  %#ok<INUSL>
data.save_det = strcmpi(mode, 'prcond');
data.save_det = true;
end
