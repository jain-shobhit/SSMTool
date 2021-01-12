function prob = coco_add_adjt(prob, fid, varargin)
%COCO_ADD_ADJT   Add adjoint to continuation problem
%
% PROB = COCO_ADD_ADJT(PROB, FID, VARARGIN)
% VARARGIN = [(@F,|@FDF,) [@DF,] DATA,] [PAR_NAMES, ['active']] OPTS...
%
% OPTS = { 'aidx' I }
% OPTS = { 'l0' L0 }
% OPTS = { 'tl0' T0 }
% OPTS = { 'adim' [N,N] }
% OPTS = 'f+df'
% OPTS = 'passchart'
% OPTS = ( 'vectorized'|'vectorised' )
% OPTS = ( 'checkderiv'|'chkdrv'|'chkderiv' [TOL] )
% OPTS = { 'remesh' @RMFUNC }
%
%   adds the adjoint associated with a previously added zero or embedded
%   monitor function to the extended continuation problem encoded in the
%   continuation problem structure PROB. If no function handle is provided,
%   then the Jacobian of the corresponding zero or monitor function is
%   used. If no function for computing the Jacobian of the adjoint is
%   provided, numerical differentiation is used. PAR_NAMES is a cell array
%   containing a list of string labels for i) the continuation parameters
%   associated with the Lagrange multipliers corresponding to the adjoints
%   of monitor functions with function type 'inactive', 'active', or
%   'internal' and ii) the continuation parameters associated with the
%   nonlinear complementarity conditions in terms of monitor functions with
%   function type 'inequality'. This list must be either empty or have
%   exactly as many names as the dimension of the vector returned by the
%   monitor function. These continuation parameters are assumed to be
%   initially inactive. This default designation can be overridden with the
%   'active' flag.
%   
% --------------------------------------------------------------
% [ [PROB] DATA [CHART] Y   ] = F  ( PROB, DATA, [CHART], U, [T] )
% [        DATA [CHART]   J ] = DF ( PROB, DATA, [CHART], U      )
% [        DATA [CHART] Y J ] = FDF( PROB, DATA, [CHART], U      )
%
% This COCO-specific function syntax applies to the evaluation of a zero or
% monitor function (F), its Jacobian (DF), or the simultaneous evaluation
% of the function and its Jacobian (FDF). 
%
%   PROB : Current continuation problem structure. Can only be returned
%     by regular or singular monitor functions. 
%
%   DATA : Current function data structure. 
%
%   CHART: Current chart.
%
%   U    : Subset of vector of continuation variables associated with
%     function dependency index set.
%
%   T    : Subset of tangent vector to current curve segment. Can only be
%     passed to regular or singular monitor functions.
%
%   Y    : Value of zero or monitor function.
%
%   J    : Value of Jacobian of zero or monitor function. Count the number
%     of input arguments in a call to FDF to determine whether this
%     variable should be evaluated.
%
% -----------------------------------------------------
% [PROB S XTR FTR] = RMFUNC(PROB, DATA, CHART, OLD_U, OLD_V)
%
% This COCO-specific function syntax applies to the adaptive remeshing of
% the adjoint of a zero or monitor function.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coco_add_adjt.m 2839 2015-03-05 17:09:01Z fschild $

%% parse input arguments

assert(isstruct(prob) && isfield(prob, 'opts') && isa(prob.opts, 'coco_opts_tree'), ...
  '%s: the first argument must be a continuation problem structure', mfilename);
assert(isfield(prob, 'efunc'), ...
  '%s: the continuation problem structure cannot be empty', mfilename);
efunc = prob.efunc;
assert(efunc.close_level<1, ...
  '%s: zero functions are closed, cannot add adjoints', mfilename);

idx   = find(strcmpi(fid, efunc.identifyers), 1);
assert(~isempty(idx), '%s: the function ''%s'' has not been defined', ...
  mfilename, fid);
[type, x0] = coco_get_func_data(prob, fid, 'type', 'x0');
assert(any(strcmp(type, {'zero' 'inactive' 'active' 'internal' ...
  'inequality'})), ...
  '%s: cannot add adjoint to function ''%s'' of type %s', ...
  mfilename, fid, type);
func  = efunc.funcs(idx);
x_idx = func.x_idx;

if ~isfield(prob, 'adjoint')
  prob.adjoint = adjoint_new([]);
  prob.complementary = complementary_new([]);
end
adjoint = prob.adjoint;

argidx = 1;

if argidx<=nargin-2 && isa(varargin{argidx}, 'function_handle')
  fhan = varargin{argidx};
  argidx = argidx + 1;
  if isa(varargin{argidx}, 'function_handle')
    dfhan = varargin{argidx};
    argidx = argidx + 1;
  else
    dfhan = [];
  end
  data   = varargin{argidx};
  argidx = argidx + 1;
else
  assert(~isempty(func.DFDX), ...
    '%s: no explicit Jacobian for corresponding zero or monitor function', mfilename);
  fhan  = func.DFDX;
  dfhan = func.DFDXDX;
  data  = func.data;
end

lnum   = 0;
lnames = {};
snum   = 0;
snames = {};
switch type
  case 'zero'
  case 'inequality'
    snames = varargin{argidx};
    if ~iscell(snames)
      snames = { snames };
    end
    if iscellstr(snames)
      argidx = argidx + 1;
      snum   = numel(snames);
      snames = reshape(snames, 1, snum);
    else
      snames = {};
    end
  otherwise
    lnames = varargin{argidx};
    if ~iscell(lnames)
      lnames = { lnames };
    end
    if iscellstr(lnames)
      argidx = argidx + 1;
      lnum   = numel(lnames);
      lnames = reshape(lnames, 1, lnum);
    else
      lnames = {};
    end
end

afopts.aidx      = [];
afopts.l0        = [];
afopts.tl0       = [];
afopts.active    = false;
afopts.vFlag     = false;
afopts.chkFlag   = false;
afopts.chkTOL    = 1e-7;
afopts.baseMode  = 1;
afopts.chartMode = 0;
afopts.remesh    = [];

while argidx<=nargin-2
  oarg   = varargin{argidx};
  argidx = argidx + 1;
  oname  = lower(oarg);
  switch oname
    
    case 'aidx'
      afopts.aidx = varargin{argidx};
      assert(isnumeric(afopts.aidx), ...
        '%s: argument %d must be an integer array', mfilename, 2+argidx);
      afopts.aidx = afopts.aidx(:)';
      argidx = argidx + 1;
      
    case 'l0'
      afopts.l0 = varargin{argidx};
      assert(isnumeric(afopts.l0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 'tl0'
      afopts.tl0 = varargin{argidx};
      assert(isnumeric(afopts.tl0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 'adim'
      afopts.adim = varargin{argidx};
      assert(isnumeric(afopts.adim) && numel(afopts.adim)==2, ...
        '%s: argument %d must be an array of two integers', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 'active'
      afopts.active = true;
      
    case 'f+df'
      afopts.baseMode = 2;
      
    case 'passchart'
      afopts.chartMode = 2;

    case { 'vectorized' 'vectorised' }
      afopts.vFlag = true;
      
    case { 'checkderiv' 'chkdrv' 'chkderiv' }
      afopts.chkFlag = true;
      if argidx<=nargin-2 && isnumeric(varargin{argidx}) ...
          && isscalar(varargin{argidx})
        afopts.chkTOL = varargin{argidx};
        argidx = argidx + 1;
      end

    case 'remesh'
      afopts.remesh = varargin{argidx};
      assert(isa(afopts.remesh, 'function_handle'), ...
        '%s: argument %d must be a function handle', mfilename, 2+argidx);
      argidx = argidx + 1;
  
    otherwise
      if ischar(oarg)
        error('%s: option ''%s'' not recognised', mfilename, oarg);
      else
        error('%s: in argument %d: expected string, got a ''%s''', ...
          mfilename, 1+argidx, class(oarg));
      end
  end
  
end

%% create function object and update adjoint

if ~isfield(afopts, 'adim')
  [data, J] = fhan(prob, data, x0);
  afopts.adim = size(J);
end
afdim = afopts.adim(1);
axdim = afopts.adim(2);
axnum = numel(afopts.aidx);
assert(axdim>=axnum, ...
  '%s: number of columns cannot be less than the size of ''aidx''', mfilename);
anum  = numel(afopts.l0);
assert(afdim==anum || anum==0, '%s: %s', mfilename, ...
  'number of initial values must equal the number of multipliers');
if anum==0
  afopts.l0  = zeros(afdim,1);
  anum       = afdim;
end
atnum = numel(afopts.tl0);
assert(atnum==0 || atnum==anum, '%s: %s', ...
  mfilename, '''tl0'' must be the same size as ''l0''');
if atnum==0
  afopts.tl0 = zeros(afdim,1);
end
if ~strcmpi(type, 'zero') && (lnum>0 || snum>0)
  assert(anum==lnum || anum==snum, '%s: %s %s', mfilename, ...
    'number of  parameter names must equal', ...
    'the number of multipliers');
end

afidx          = adjoint.a_dim(1) + (1:afdim);
axidx          = adjoint.a_dim(2) + (1:axdim-axnum);
adjoint.a_dim  = adjoint.a_dim + [afdim, axdim-axnum];
adjoint.af_idx = [ adjoint.af_idx afidx ];
adjoint.ax_idx = [ adjoint.ax_idx axidx ];
axidx          = [ afopts.aidx    axidx ];

func.identifyer = fid;
func.F          = fhan;
func.DFDX       = dfhan;
func.data       = data;
func.lnum       = lnum;
func.lnames     = lnames;
func.snum       = snum;
func.snames     = snames;
func.af_idx     = afidx;
func.ax_idx     = axidx;
func.x_idx      = x_idx;
func.vectorised = afopts.vFlag;
func.chkdrv     = afopts.chkFlag;
func.chkTOL     = afopts.chkTOL;
func.call_mode  = afopts.baseMode + afopts.chartMode;
func.remesh     = afopts.remesh;

adjoint.funcs       = [ adjoint.funcs func ];
adjoint.identifyers = [ adjoint.identifyers fid ];

if idx > numel(adjoint.midx)
  adjoint.midx{idx} = numel(adjoint.identifyers);
else
  adjoint.midx{end} = [adjoint.midx{end} numel(adjoint.identifyers)];
end

adjoint.l0  = [ adjoint.l0  ;  afopts.l0 ];
adjoint.tl0 = [ adjoint.tl0 ; afopts.tl0 ];
switch type
  case 'zero'
  case 'inequality'
    if snum>0
      if afopts.active
        adjoint.sactive   = [ adjoint.sactive   ; numel(adjoint.s_idx)+(1:afdim)' ];
      else
        adjoint.sinactive = [ adjoint.sinactive ; numel(adjoint.s_idx)+(1:afdim)' ];
      end
      adjoint.snames = [ adjoint.snames        snames ];
    end
    adjoint.s_idx  = [ adjoint.s_idx ; func.af_idx' ];
  otherwise
    if lnum>0
      if afopts.active
        adjoint.lactive   = [ adjoint.lactive   ; numel(adjoint.l_idx)+(1:afdim)' ];
      else
        adjoint.linactive = [ adjoint.linactive ; numel(adjoint.l_idx)+(1:afdim)' ];
      end
      adjoint.lnames = [ adjoint.lnames        lnames ];
    end
    adjoint.l_idx  = [ adjoint.l_idx ; func.af_idx' ];
end

prob.adjoint = adjoint;

end
