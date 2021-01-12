function prob = coco_add_func(prob, fid, varargin)
%COCO_ADD_FUNC   Basic COCO constructor - add zero or monitor function
%
% PROB = COCO_ADD_FUNC(PROB, FID, VARARGIN)
% VARARGIN = (@F, [@DF, [@DDF,]] | @FDF, [@DDF,]) DATA, TYPE_SPEC, OPTS...
%
% TYPE_SPEC = { 'zero' }
% TYPE_SPEC = { ('inactive'|'active'|'internal') PAR_NAMES }
% TYPE_SPEC = { ('regular'|'singular') PAR_NAMES }
% TYPE_SPEC = { 'inequality' PAR_NAMES }
%
% OPTS = { ('uidx'|'xidx') ('all'|I) }
% OPTS = { ('u0'|'x0') U0 }
% OPTS = { 't0' T }
% OPTS = { 'j0' J }
% OPTS = { ('f0'|'y0') F }
% OPTS = { ('fdim'|'ydim') N }
% OPTS = 'f+df'
% OPTS = 'passchart'
% OPTS = 'passtangent'
% OPTS = 'returnsprob'
% OPTS = ( 'vectorized'|'vectorised' )
% OPTS = ( 'checkderiv'|'chkdrv'|'chkderiv' [TOL] )
% OPTS = { 'remesh' @RMFUNC }
% OPTS = { 'copy' @CPFUNC }
% OPTS = { 'requires' FID_LIST }
% OPTS = { ('prop'|'property') PROP_NAME PROP_VAL }
%
%   adds a new function to the extended continuation problem encoded in the
%   continuation problem structure PROB. If no function for computing the
%   Jacobian is present, numerical differentiation is used. TYPE defines
%   how this function is added to the system; see explanations below.
%   PAR_NAMES is a cell array containing a list of string labels for the
%   continuation parameters associated to a monitor function. This list
%   must have exactly as many names as the dimension of the vector returned
%   by the monitor function.
%
%   Valid values for TYPE are:
%
%   'zero' : add the function to the zero problem of the extended
%      continuation problem.
%
%   'inactive' : add the function as an embedded monitor function to the
%      extended continuation problem and mark the associated continuation
%      parameters as inactive. Inactive continuation parameters impose
%      additional constraints during continuation. Initially inactive
%      continuation parameters can be activated by an exchange with an
%      active continuation parameter using coco_xchg_pars or by explicitly
%      releasing them in the call to the coco entry-point function.
%
%   'active' : add the function as an embedded monitor function to the
%      extended continuation problem and mark the associated continuation
%      parameters as active. Active continuation parameters track the
%      values of monitor functions during continuation and are checked for
%      events. Active parameters can be inactivated by an exchange with an
%      inactive continuation parameter using coco_xchg_pars.
%
%   'internal' : add the function as an embedded monitor function to the
%      extended continuation problem and mark the associated continuation
%      parameters as internal. Internal continuation parameters are like
%      active continuation parameters with the additional property that
%      they will be exchanged automatically if the user overspecifies
%      continuation parameters, that is, specifies more initially inactive
%      continuation parameters than the dimensional deficit. Internal
%      parameters are automatically included with the screen output during
%      continuation.
%
%   'inequality' : add the function as an embedded monitor function to the
%      extended continuation problem and mark the associated continuation
%      parameter as active. The corresponding Lagrange multipliers are
%      treated as if the monitor function were a zero function. Each
%      inequality monitor function and corresponding Lagrange multiplier
%      contribute to the nonlinear complementarity constraints in terms of
%      a "relaxation" continuation parameter that is nonzero when the
%      corresponding complementarity condition is violated.
%
%   'regular' : add the function to the set of monitor functions that are
%      not embedded in the continuation problem and treat events associated
%      to the corresponding continuation parameters as regular special
%      points (IFT is satisfied). Regular continuation parameters do not
%      increase the dimensional deficit of a continuation problem.
%
%   'singular' : add the function to the set of monitor functions that are
%      not embedded in the continuation problem and treat events associated
%      to the corresponding continuation parameters as singular special
%      points (IFT may be violated). Singular continuation parameters do
%      not increase the dimensional deficit of a continuation problem.
%
%   The 'all' option associated with the 'uidx'|'xidx' flag allows for
%   the definition of a zero or monitor function that depends on the entire
%   vector of unknown variables (or the corresponding tangent vector) in a
%   continuation problem. The definition of functions associated with the
%   'all' option must be accompanied by values for either the 'fdim' flag
%   or the 'f0'|'y0' flag, cannot be associated with an explicit remesh
%   function, and cannot introduce new variables.
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
%   CHART : Current chart.
%
%   U : Subset of vector of continuation variables associated with function
%     dependency index set.
%
%   T : Subset of tangent vector to current curve segment associated with
%     function dependency index set. Can only be passed to regular or
%     singular monitor functions.
%
%   Y : Value of zero or monitor function.
%
%   J : Value of Jacobian of zero or monitor function. Count the number
%     of output arguments in a call to FDF to determine whether this
%     variable should be evaluated.
%
% -----------------------------------------------------
% [PROB S TR] = RMFUNC(PROB, DATA, CHART, OLD_U, OLD_V)
%
% This COCO-specific function syntax applies to the adaptive remeshing of a
% zero or monitor function.
%
%   PROB : Current continuation problem structure. 
%
%   DATA : Current function data structure. 
%
%   CHART : Current chart.
%
%   OLD_U : Subset of vector of continuation variables associated with
%     function dependency index set prior to adaptive remeshing.
%
%   OLD_V : Subset of matrix of tangent vectors associated with function
%     dependency index set prior to adaptive remeshing.
%
%   S : String output indicating successful remeshing ('success'), failed
%     remeshing ('failure') or the need for further remeshing ('repeat').
%
%   TR : Integer vector of length equal to size of set difference of
%     function dependency index set and union of function dependency index
%     sets of previously defined functions prior to adaptive remeshing.
%     Nonzero entries indexed by mesh-invariant elements and equal to
%     corresponding indices after adaptive remeshing.
%
% -------------------------
% DATA = CPFUNC(PROB, DATA)
%
% This COCO-specific function syntax applies to the copying of function
% data.
%
%   PROB : Current continuation problem structure.
%
%   DATA : Current function data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coco_add_func.m 3099 2019-06-07 17:31:05Z hdankowicz $

%% parse input arguments

assert(isstruct(prob) && isfield(prob, 'opts') && isa(prob.opts, 'coco_opts_tree'), ...
  '%s: the first argument must be a continuation problem structure', mfilename);
if ~isfield(prob, 'efunc') 
  prob.efunc = efunc_new([]);
end
efunc = prob.efunc;

coco_opts_tree.check_path(fid);
check_names(['efunc' 'mfunc' efunc.identifyers], fid, 'function');

argidx = 1;

fhan = varargin{argidx};
assert(isa(fhan, 'function_handle'), ...
  '%s: argument %d must be a function handle', mfilename, 2+argidx);
argidx = argidx + 1;

if isa(varargin{argidx}, 'function_handle')
  dfhan = varargin{argidx};
  argidx = argidx + 1;
else
  dfhan = [];
end

if isa(varargin{argidx}, 'function_handle')
  ddfhan = varargin{argidx};
  argidx = argidx + 1;
else
  ddfhan = [];
end

data   = varargin{argidx};
type   = lower(varargin{argidx+1});
argidx = argidx + 2;

switch type
  case 'zero'
    assert(efunc.close_level<1, ...
      '%s: zero functions are closed, cannot add function ''%s''', ...
      mfilename, fid);
  case { 'inactive' 'active' 'internal' 'inequality' }
    assert(efunc.close_level<1, ...
      '%s: embedded monitor functions are closed, cannot add function ''%s''', ...
      mfilename, fid);
  case { 'regular' 'singular' }
    assert(efunc.close_level<2, ...
      '%s: monitor functions are closed, cannot add function ''%s''', ...
      mfilename, fid);
  otherwise
    error('%s: type ''%s'' not recognised', mfilename, type);
end

switch type
  case 'zero'
    funcarray = 'zero';
  case { 'inactive' 'active' 'internal' 'inequality' }
    funcarray = 'embedded';
  case 'regular'
    funcarray = 'regular';
  case 'singular'
    funcarray = 'singular';
end

switch type % associate function with continuation parameters
  case 'zero'
    pnum   = 0;
    pnames = {};
  otherwise
    pnames = varargin{argidx};
    if ~iscell(pnames)
      pnames = { pnames };
    end
    assert(iscellstr(pnames), ...
      '%s: argument %d must be a cell array of strings', mfilename, 2+argidx);
    argidx = argidx + 1;
    pnum   = numel(pnames);
    pnames = reshape(pnames, 1, pnum);
end

efopts.vFlag     = false;
efopts.chkFlag   = false;
efopts.chkTOL    = 1e-7;
efopts.baseMode  = 1;
efopts.chartMode = 0;
efopts.tanMode   = 0;
efopts.optsMode  = 0;
efopts.x0        = [];
efopts.xidx      = [];
efopts.t0        = [];
efopts.remesh    = [];
efopts.copy      = [];
efopts.requires  = {};
efopts.props     = struct();

while argidx<=nargin-2
  oarg   = varargin{argidx};
  argidx = argidx + 1;
  oname  = lower(oarg);
  switch oname
    
    case { 'uidx' 'xidx' }
      efopts.xidx = varargin{argidx};
      assert(isnumeric(efopts.xidx) || strcmpi(efopts.xidx, 'all'), ...
        '%s: argument %d must be an integer array', mfilename, 2+argidx);
      efopts.xidx = efopts.xidx(:)';
      argidx = argidx + 1;
      
    case { 'u0' 'x0' }
      assert(any(strcmp(type, ...
        {'zero' 'inactive' 'active' 'internal' 'inequality'})), ...
        '%s: %s: %s functions cannot add continuation variables', ...
        mfilename, fid, type);
      efopts.x0 = varargin{argidx};
      assert(isnumeric(efopts.x0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 't0'
      assert(any(strcmp(type, ...
        {'zero' 'inactive' 'active' 'internal' 'inequality'})), ...
        '%s: %s: %s functions cannot define tangent speeds', ...
        mfilename, fid, type);
      efopts.t0 = varargin{argidx};
      assert(isnumeric(efopts.t0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 'j0'
      assert(any(strcmp(type, ...
        {'inactive' 'active' 'internal' 'inequality'})), ...
        '%s: %s: %s functions cannot define differential', ...
        mfilename, fid, type);
      efopts.J0 = varargin{argidx};
      assert(isnumeric(efopts.J0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case { 'f0' 'y0' }
      efopts.fdim = numel(varargin{argidx});
      argidx = argidx + 1;
      
    case { 'fdim' 'ydim' }
      efopts.fdim = varargin{argidx};
      assert(isnumeric(efopts.fdim) && isscalar(efopts.fdim), ...
        '%s: argument %d must be an integer', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 'f+df'
      efopts.baseMode = 2;
      
    case 'passchart'
      efopts.chartMode = 2;
      
    case 'passtangent'
      assert(any(strcmp(type, {'regular' 'singular'})), ...
        '%s: error when adding function ''%s'',\n%s', mfilename, fid, ...
        'tangent can only be passed to regular or singular monitor functions');
      efopts.tanMode = 4;
      
    case { 'returnsprob' 'returnsopts' }
      assert(any(strcmp(type, {'regular' 'singular'})), ...
        '%s: error when adding function ''%s'',\n%s %s', mfilename, fid, ...
        'continuation problem structure can only be returned by', ...
        'regular or singular monitor functions');
      efopts.optsMode = 8;

    case { 'vectorized' 'vectorised' }
      efopts.vFlag = true;
      
    case {'checkderiv' 'chkdrv' 'chkderiv' }
      efopts.chkFlag = true;
      if argidx<=nargin-2 && isnumeric(varargin{argidx}) ...
          && isscalar(varargin{argidx})
        efopts.chkTOL = varargin{argidx};
        argidx = argidx + 1;
      end

    case 'remesh'
      efopts.remesh = varargin{argidx};
      assert(isa(efopts.remesh, 'function_handle'), ...
        '%s: argument %d must be a function handle', mfilename, 2+argidx);
      argidx = argidx + 1;
      
    case 'copy'
      efopts.copy = varargin{argidx};
      assert(isa(efopts.copy, 'function_handle'), ...
        '%s: argument %d must be a function handle', mfilename, 2+argidx);
      argidx = argidx + 1;

    case 'requires'
      efopts.requires = varargin{argidx};
      if ~iscell(efopts.requires)
        efopts.requires = { efopts.requires };
      end
      idx = cellfun(@(x)any(strcmpi(x,efunc.identifyers)),efopts.requires);
      assert(all(idx), '%s: required function ''%s'' not found', ...
        mfilename, efopts.requires{find(~idx,1)});
      idx = cellfun(@(x)find(strcmpi(x,efunc.identifyers)),efopts.requires);
      switch type
        case { 'inactive' 'active' 'internal' 'inequality'}
          flag = all(ismember(idx, efunc.embedded));
          assert(flag, '%s: %s %s', mfilename, ...
            'embedded monitor functions can only require other', ...
            'embedded monitor functions');
        case 'zero'
          flag = all(ismember(idx, [efunc.embedded, efunc.zero]));
          assert(flag, '%s: %s %s', mfilename, ...
            'zero functions can only require other zero or', ...
            'embedded monitor functions');
        case 'regular'
          flag = all(ismember(idx, [efunc.embedded, efunc.regular]));
          assert(flag, '%s: %s %s', mfilename, ...
            'regular monitor functions can only require other', ...
            'regular or embedded monitor functions');
        case 'singular'
          flag = all(ismember(idx, [efunc.embedded, efunc.regular efunc.singular]));
          assert(flag, '%s: %s %s', mfilename, ...
            'singular monitor functions can only require other', ...
            'singular, regular, or embedded monitor functions');
      end
      argidx = argidx + 1;
      
    case { 'prop' 'property' }
      name   = varargin{argidx};
      val    = varargin{argidx+1};
      argidx = argidx + 2;
      efopts.props.(name) = val;
      
    otherwise
      if ischar(oarg)
        error('%s: option ''%s'' not recognised', mfilename, oarg);
      else
        error('%s: in argument %d: expected string, got a ''%s''', ...
          mfilename, 1+argidx, class(oarg));
      end
  end
  
end

xnum  = numel(efopts.x0);
if strcmpi(efopts.xidx, 'all')
  assert(xnum==0, '%s: %s %s', mfilename, 'the option ''u0|x0''', ...
    'cannot be included when ''uidx|xidx'' equals ''all''');
  assert(isfield(efopts, 'fdim') && efopts.fdim>0, '%s: %s %s', ...
    mfilename, 'either of the options ''f0|y0'' or ''fdim|ydim''', ...
    'must be included when ''uidx|xidx'' equals ''all''');
  assert(~isa(efopts.remesh, 'function_handle'), '%s: %s %s', ...
    mfilename, 'a remesh function cannot be included', ...
    'when ''uidx|xidx'' equals ''all''');
else
  if any(strcmp(type, {'zero' 'internal' 'inactive' 'active' 'inequality'}))
    assert(xnum+numel(efopts.xidx)>0, '%s: %s %s', mfilename, ...
      'at least one of ''uidx|xidx'' or ''u0|x0''', ...
      'must be present and non-empty');
  end
  assert(isempty(efopts.t0) || numel(efopts.t0)==xnum, '%s: %s', ...
    mfilename, '''t0'' must either be empty or the same size as ''x0''');
  x0          = efopts.x0(:);
  if ~isempty(efopts.t0)
    t0        = efopts.t0(:);
  else
    t0        = zeros(xnum,1);
  end
  xidx        = efunc.x_dim + (1:xnum);
  efunc.x_dim = efunc.x_dim + xnum;
  efunc.x0    = [ efunc.x0              ;   x0 ];
  x0          = [ efunc.x0(efopts.xidx) ;   x0 ];
  efunc.tx    = [ efunc.tx              ;   t0 ];
  efunc.x_idx = [ efunc.x_idx             xidx ];
  efopts.xidx = [ efopts.xidx             xidx ];
end

%% create function object and update efunc.f_dim

func.identifyer   = fid;
func.F            = fhan;
func.DFDX         = dfhan;
func.DFDXDX       = ddfhan;
func.data         = data;
func.type         = type;
func.pnum         = pnum;
func.pnames       = pnames;
func.x_idx        = efopts.xidx;
func.vectorised   = efopts.vFlag;
func.chkdrv       = efopts.chkFlag;
func.chkTOL       = efopts.chkTOL;
func.call_mode    = efopts.baseMode + efopts.chartMode ...
  + efopts.tanMode + efopts.optsMode;
func.remesh       = efopts.remesh;
func.copy         = efopts.copy;
func.requires     = efopts.requires;
func.props        = efopts.props;

switch type
  case {'zero' 'regular' 'singular'}
    if ~isfield(efopts, 'fdim') % don't call F if f0 or fdim given by user
      [prob, func.data, efunc.chart, f] = ...
        efunc_call_F(prob, func.data, efunc.chart, func, x0);
      efopts.fdim = numel(f);
    end
  case {'inactive' 'active' 'internal' 'inequality'}
    if isfield(efopts, 'J0')
      J = efopts.J0;
      if ~isfield(efopts, 'fdim') % don't call F if f0 or fdim given by user
        [prob, func.data, efunc.chart, f] = ...
          efunc_call_F(prob, func.data, efunc.chart, func, x0);
        efopts.fdim = numel(f);
      end
    elseif isfield(efopts, 'fdim')
      [prob, func.data, efunc.chart, J] = ...
        efunc_call_DF(prob, func.data, efunc.chart, func, x0);
    else
      [prob, func.data, efunc.chart, f, J] = ...
        efunc_call_FDF(prob, func.data, efunc.chart, func, x0);
      efopts.fdim = numel(f);
    end
end

fdim = efopts.fdim;
switch type
  case { 'zero' }
    fidx             = efunc.f_dim + (1:fdim);
    efunc.f_dim      = efunc.f_dim + fdim;
    efunc.pidx2fidx  = [ efunc.pidx2fidx ; fidx' ];
    efunc.tp(fidx,1) = 0;
  case { 'inactive' 'active' 'internal' 'inequality'}
    fidx             = efunc.f_dim + (1:fdim);
    efunc.f_dim      = efunc.f_dim + fdim;
    efunc.pidx2fidx  = [ efunc.pidx2fidx ; fidx' ];
    efunc.tp(fidx,1) = J*efunc.tx(efopts.xidx);
  otherwise
    fidx            = [];
    efunc.pidx2fidx = [ efunc.pidx2fidx ; zeros(fdim,1) ];
end
func.f_idx  = fidx;
efunc.f_idx = [efunc.f_idx fidx];

switch type
  case 'zero'
    midx            = [];
    pnames          = cell(1,fdim);
    efunc.pidx2midx = [ efunc.pidx2midx ; zeros(fdim,1) ];
  otherwise
    if fdim~=pnum
      error('%s: number of parameter names must match dimension of function', ...
        mfilename);
    end
    midx            = efunc.m_dim + (1:fdim);
    efunc.m_dim     = efunc.m_dim + fdim;
    efunc.pidx2midx = [ efunc.pidx2midx ; midx' ];
end
func.m_idx          = midx;

func.req_idx = midx;
for i=1:numel(func.requires)
  idx = strcmpi(func.requires{i}, efunc.identifyers);
  efunc.funcs(idx).req_idx = [efunc.funcs(idx).req_idx midx];
end

efunc.funcs       = [ efunc.funcs func ];
efunc.(funcarray) = [ efunc.(funcarray) numel(efunc.funcs) ];
efunc.identifyers = [ efunc.identifyers fid ];

%% add pnames to list of parameters and update efunc.p_dim
%  check for duplicates, add names and compute parameter indices

if pnum>0
  check_names(efunc.idx2par, pnames, 'parameter');
end
idx              = numel(efunc.idx2par) + (1:pnum);
efunc.idx2par    = [ efunc.idx2par pnames ];
pararray         = sprintf('%s_pars', type);
efunc.(pararray) = [ efunc.(pararray) idx ];

switch type
  case { 'active' 'internal' 'inequality' }
    efunc.p_dim = efunc.p_dim + pnum;
end

%% return efunc in prob

prob.efunc = efunc;

%% add pending functions

if prob.efunc.add_pending
  prob = coco_add_pending(prob);
end

end

%%
function check_names(list, names, what)
%CHECK_NAMES   check for duplicate names in a list
%
%   CHECK_NAMES(LIST, NAMES, WHAT) throws an error if one of the strings in
%   cell array names is already present in cell array list. The function
%   returns if no duplicates are found. The error message will contain the
%   string WHAT for readability.

if ~iscell(names)
  names = { names };
end

for i=1:numel(names)
  name = names{i};
  if isempty(name)
    error('%s: %s names must not be empty', mfilename, what);
  end
  if any(strcmpi(name, list))
    error('%s: %s with name ''%s'' already defined', mfilename, what, name);
  end
end

end
