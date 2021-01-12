function prob = coco_add_comp(prob, fid, varargin)
%COCO_ADD_COMP   Basic COCO constructor - add complementary zero or monitor function
%
% PROB = COCO_ADD_COMP(PROB, FID, VARARGIN)
% VARARGIN = (@F, [@DF,] | @FDF,) DATA, TYPE_SPEC, OPTS...
%
% TYPE_SPEC = { 'zero' }
% TYPE_SPEC = { ('inactive'|'active'|'internal') PAR_NAMES }
%
% OPTS = { 'uidx' I }
% OPTS = { 'lidx' I }
% OPTS = { 'vidx' I }
% OPTS = { 'v0'  V }
% OPTS = { 'tv0' T }
% OPTS = { 'j0'  J }
% OPTS = { 'f0'  F }
% OPTS = { 'fdim' N }
% OPTS = 'f+df'
% OPTS = ( 'vectorized'|'vectorised' )
% OPTS = ( 'checkderiv'|'chkdrv'|'chkderiv' [TOL] )
% OPTS = { 'remesh' @RMFUNC }
% OPTS = { 'copy' @CPFUNC }
% OPTS = { 'requires' FID_LIST }
% OPTS = { ('prop'|'property') PROP_NAME PROP_VAL }
%
%   adds a new complementary function to the augmented continuation problem
%   encoded in the continuation problem structure PROB. If no function for
%   computing the Jacobian is present, numerical differentiation is used.
%   TYPE defines how this function is added to the system; see explanations
%   below. PAR_NAMES is a cell array containing a list of string labels for
%   the continuation parameters associated to a monitor function. This list
%   must have exactly as many names as the dimension of the vector returned
%   by the monitor function.
%
%   Valid values for TYPE are:
%
%   'zero' : add the function to the complementary zero problem of the
%      continuation problem.
%
%   'inactive' : add the function as an embedded complementary monitor
%      function to the continuation problem and mark the associated
%      complementary continuation parameters as inactive. Inactive
%      continuation parameters impose additional constraints during
%      continuation. Initially inactive continuation parameters can be
%      activated by an exchange with an active continuation parameter using
%      coco_xchg_pars or by explicitly releasing them in the call to the
%      coco entry-point function.
%
%   'active' : add the function as an embedded complementary monitor
%      function to the continuation problem and mark the associated
%      complementary continuation parameters as active. Active continuation
%      parameters track the values of monitor functions during continuation
%      and are checked for events. Active parameters can be inactivated by
%      an exchange with an inactive continuation parameter using
%      coco_xchg_pars.
%
%   'internal' : add the function as an embedded complementary monitor
%      function to the extended continuation problem and mark the
%      associated complementary continuation parameters as internal.
%      Internal continuation parameters are like active continuation
%      parameters with the additional property that they will be exchanged
%      automatically if the user overspecifies continuation parameters,
%      that is, specifies more initially inactive continuation parameters
%      than the dimensional deficit. Internal parameters are automatically
%      included with the screen output during continuation.
%
%   'regular' : add the function to the set of complementary monitor
%      functions that are not embedded in the continuation problem and
%      treat events associated to the corresponding continuation parameters
%      as regular special points (IFT is satisfied). Regular continuation
%      parameters do not increase the dimensional deficit of a continuation
%      problem.
%
%   'singular' : add the function to the set of complementary monitor
%      functions that are not embedded in the continuation problem and
%      treat events associated to the corresponding continuation parameters
%      as singular special points (IFT may be violated). Singular
%      continuation parameters do not increase the dimensional deficit of a
%      continuation problem.
%
% --------------------------------------------------------------
% [ DATA Y   ] = F  ( PROB, DATA, U, L, V )
% [ DATA   J ] = DF ( PROB, DATA, U, L, V )
% [ DATA Y J ] = FDF( PROB, DATA, U, L, V )
%
% This COCO-specific function syntax applies to the evaluation of a zero or
% monitor function (F), its Jacobian (DF), or the simultaneous evaluation
% of the function and its Jacobian (FDF). 
%
%   PROB : Current continuation problem structure.
%
%   DATA : Current function data structure. 
%
%   U : Subset of vector of continuation variables associated with function
%     dependency index set.
%
%   L : Subset of vector of continuation multipliers associated with
%     function dependency index set.
%
%   V : Subset of vector of complementary continuation variables associated
%     with function dependency index set.
%
%   Y : Value of complementary zero or monitor function.
%
%   J : Value of Jacobian of complementary zero or monitor function. Count
%     the number of output arguments in a call to FDF to determine whether
%     this variable should be evaluated.
%
% -----------------------------------------------------
% [PROB S TR] = RMFUNC(PROB, DATA, CHART, OLD_U, OLD_V)
%
% This COCO-specific function syntax applies to the adaptive remeshing of a
% complementary zero or monitor function.
%
%   PROB : Current continuation problem structure. 
%
%   DATA : Current function data structure. 
%
%   CHART : Current chart.
%
%   OLD_U : Subset of vector of complementary continuation variables
%     associated with function dependency index set prior to adaptive
%     remeshing.
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

% Copyright (C) Harry Dankowicz, Mingwu Li

%% parse input arguments

assert(isstruct(prob) && isfield(prob, 'opts') && isa(prob.opts, 'coco_opts_tree'), ...
  '%s: the first argument must be a continuation problem structure', mfilename);
assert(isfield(prob, 'efunc'), ...
  '%s: the continuation problem structure cannot be empty', mfilename);
efunc = prob.efunc;
assert(isfield(prob, 'adjoint'), ...
  '%s: the continuation problem structure must contain adjoint functions', mfilename);
adjoint = prob.adjoint;

if ~isfield(prob, 'complementary')
  prob.complementary = complementary_new([]);
end
complementary = prob.complementary;

coco_opts_tree.check_path(fid);
check_names(['complementary' 'efunc' 'mfunc' complementary.identifyers ...
  efunc.identifyers], fid, 'function');

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

data   = varargin{argidx};
type   = lower(varargin{argidx+1});
argidx = argidx + 2;

switch type
  case 'zero'
    assert(efunc.close_level<1, ...
      '%s: zero functions are closed, cannot add complementary function ''%s''', ...
      mfilename, fid);
  case { 'inactive' 'active' 'internal'}
    assert(efunc.close_level<1, ...
      '%s: embedded monitor functions are closed, cannot add complementary function ''%s''', ...
      mfilename, fid);
  case { 'regular' 'singular' }
    assert(efunc.close_level<2, ...
      '%s: monitor functions are closed, cannot add complementary function ''%s''', ...
      mfilename, fid);
  otherwise
    error('%s: type ''%s'' not recognised', mfilename, type);
end

switch type
  case 'zero'
    funcarray = 'zero';
  case { 'inactive' 'active' 'internal' }
    funcarray = 'embedded';
  case 'regular'
    funcarray = 'regular';
  case 'singular'
    funcarray = 'singular';
end

args = { type };
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
    args = [ args, {pnames} ];
end

cfopts.vFlag    = false;
cfopts.chkFlag  = false;
cfopts.chkTOL   = 1e-7;
cfopts.baseMode = 1;
cfopts.uidx     = [];
cfopts.lidx     = [];
cfopts.vidx     = [];
cfopts.v0       = [];
cfopts.tv0      = [];
cfopts.remesh   = [];
cfopts.copy     = [];
cfopts.requires = {};
cfopts.props    = struct();

while argidx<=nargin-2
  oarg   = varargin{argidx};
  argidx = argidx + 1;
  oname  = lower(oarg);
  switch oname
    
    case 'uidx'
      cfopts.uidx = varargin{argidx};
      assert(isnumeric(cfopts.uidx), ...
        '%s: argument %d must be an integer array', mfilename, 2+argidx);
      cfopts.uidx = cfopts.uidx(:)';
      argidx = argidx + 1;
      
    case 'lidx'
      cfopts.lidx = varargin{argidx};
      assert(isnumeric(cfopts.lidx), ...
        '%s: argument %d must be an integer array', mfilename, 2+argidx);
      cfopts.lidx = cfopts.lidx(:)';
      argidx = argidx + 1;
      
    case 'vidx'
      cfopts.vidx = varargin{argidx};
      assert(isnumeric(cfopts.vidx), ...
        '%s: argument %d must be an integer array', mfilename, 2+argidx);
      cfopts.vidx = cfopts.vidx(:)';
      argidx = argidx + 1;
      
    case 'v0'
      cfopts.v0 = varargin{argidx};
      assert(isnumeric(cfopts.v0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      args   = [ args, {'u0', cfopts.v0} ]; %#ok<AGROW>
      
    case 'tv0'
      cfopts.tv0 = varargin{argidx};
      assert(isnumeric(cfopts.tv0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      args   = [ args, {'t0', cfopts.tv0} ]; %#ok<AGROW>
      
    case 'j0'
      assert(any(strcmp(type, ...
        {'inactive' 'active' 'internal'})), ...
        '%s: %s: %s functions cannot define differential', ...
        mfilename, fid, type);
      cfopts.J0 = varargin{argidx};
      assert(isnumeric(cfopts.J0), ...
        '%s: argument %d must be a numeric array', mfilename, 2+argidx);
      argidx = argidx + 1;
      args   = [ args, {'j0', cfopts.J0} ]; %#ok<AGROW>
      
    case 'f0'
      cfopts.fdim = numel(varargin{argidx});
      argidx = argidx + 1;
      args   = [ args, {'fdim', cfopts.fdim} ]; %#ok<AGROW>
      
    case 'fdim'
      cfopts.fdim = varargin{argidx};
      assert(isnumeric(cfopts.fdim) && isscalar(cfopts.fdim), ...
        '%s: argument %d must be an integer', mfilename, 2+argidx);
      argidx = argidx + 1;
      args   = [ args, {'fdim', cfopts.fdim} ]; %#ok<AGROW>
      
    case 'f+df'
      cfopts.baseMode = 2;

    case 'vectorized'
      cfopts.vFlag = true;
      args   = [ args, 'vectorized' ]; %#ok<AGROW>
      
    case {'checkderiv' 'chkdrv' 'chkderiv' }
      cfopts.chkFlag = true;
      if argidx<=nargin-2 && isnumeric(varargin{argidx}) ...
          && isscalar(varargin{argidx})
        cfopts.chkTOL = varargin{argidx};
        argidx = argidx + 1;
      end
      args   = [ args, {'checkderiv', cfopts.chkTOL} ]; %#ok<AGROW>

    case 'remesh'
      cfopts.remesh = varargin{argidx};
      assert(isa(cfopts.remesh, 'function_handle'), ...
        '%s: argument %d must be a function handle', mfilename, 2+argidx);
      argidx = argidx + 1;
      args   = [ args, {'remesh', cfopts.remesh} ]; %#ok<AGROW>
      
    case 'copy'
      cfopts.copy = varargin{argidx};
      assert(isa(cfopts.copy, 'function_handle'), ...
        '%s: argument %d must be a function handle', mfilename, 2+argidx);
      argidx = argidx + 1;
      args   = [ args, {'copy', cfopts.copy} ]; %#ok<AGROW>

    case 'requires'
      cfopts.requires = varargin{argidx};
      if ~iscell(cfopts.requires)
        cfopts.requires = { cfopts.requires };
      end
      idx = cellfun(@(x) any(strcmpi(x, ...
        [efunc.identifyers complementary.identifyers])),cfopts.requires);
      assert(all(idx), '%s: required function ''%s'' not found', ...
        mfilename, cfopts.requires{find(~idx,1)});
      idx = cellfun(@(x) find(strcmpi(x,...
        [efunc.identifyers complementary.identifyers])),cfopts.requires);
      switch type
        case { 'inactive' 'active' 'internal' 'inequality'}
          flag = all(ismember(idx, [efunc.embedded, complementary.embedded]));
          assert(flag, '%s: %s %s', mfilename, ...
            'embedded monitor functions can only require other', ...
            'embedded monitor functions');
        case 'zero'
          flag = all(ismember(idx, [efunc.embedded, efunc.zero, ...
            complementary.embedded, complementary.zero]));
          assert(flag, '%s: %s %s', mfilename, ...
            'zero functions can only require other zero or', ...
            'embedded monitor functions');
        case 'regular'
          flag = all(ismember(idx, [efunc.embedded, efunc.regular, ...
            complementary.embedded, complementary.regular]));
          assert(flag, '%s: %s %s', mfilename, ...
            'regular monitor functions can only require other', ...
            'regular or embedded monitor functions');
        case 'singular'
          flag = all(ismember(idx, [efunc.embedded, efunc.regular, ...
            efunc.singular, complementary.embedded, complementary.regular, ...
            complementary.singular]));
          assert(flag, '%s: %s %s', mfilename, ...
            'singular monitor functions can only require other', ...
            'singular, regular, or embedded monitor functions');
      end
      argidx = argidx + 1;
      args   = [ args, {'requires', cfopts.requires} ]; %#ok<AGROW>
      
    case { 'prop' 'property' }
      name   = varargin{argidx};
      val    = varargin{argidx+1};
      argidx = argidx + 2;
      cfopts.props.(name) = val;
      args   = [ args, {'prop', cfopts.props} ]; %#ok<AGROW>
      
    otherwise
      if ischar(oarg)
        error('%s: option ''%s'' not recognised', mfilename, oarg);
      else
        error('%s: in argument %d: expected string, got a ''%s''', ...
          mfilename, 1+argidx, class(oarg));
      end
  end
  
end

u0 = efunc.x0(cfopts.uidx);
l0 = adjoint.l0(cfopts.lidx);

vnum = numel(cfopts.v0);
assert(numel(u0)+numel(l0)+vnum+numel(cfopts.vidx)>0, ...
  '%s: %s %s', mfilename, 'at least one of ''uidx'', ''lidx'',', ...
  '''vidx'' or ''v0'' must be present and non-empty');
assert(isempty(cfopts.tv0) || numel(cfopts.tv0)==vnum, '%s: %s', ...
  mfilename, '''tv0'' must either be empty or the same size as ''v0''');
v0    = cfopts.v0(:);
if ~isempty(cfopts.tv0)
  tv0 = cfopts.tv0(:);
else
  tv0 = zeros(vnum,1);
end

if vnum>0
  vidx = complementary.v_dim + (1:vnum);
else
  vidx = [];
end
complementary.v_dim = complementary.v_dim + vnum;
complementary.v0    = [ complementary.v0              ;   v0 ];
v0                  = [ complementary.v0(cfopts.vidx) ;   v0 ];
complementary.tv    = [ complementary.tv              ;  tv0 ];
complementary.v_idx = [ complementary.v_idx             vidx ];
cfopts.vidx         = [ cfopts.vidx                     vidx ];

%% create function object and update complementary.f_dim

func.identifyer = fid;
func.F          = fhan;
func.DFDX       = dfhan;
func.type       = type;
func.pnum       = pnum;
func.pnames     = pnames;
func.u_idx      = cfopts.uidx;
func.l_idx      = cfopts.lidx;
func.v_idx      = cfopts.vidx;
func.v_num      = vnum;
func.vectorised = cfopts.vFlag;
func.chkdrv     = cfopts.chkFlag;
func.chkTOL     = cfopts.chkTOL;
func.call_mode  = cfopts.baseMode;
func.remesh     = cfopts.remesh;
func.requires   = cfopts.requires;
func.copy       = cfopts.copy;
func.props      = cfopts.props;

if ~isfield(cfopts, 'fdim') % don't call F if f0 or fdim given by user
  [data, f] = complementary_call_F(prob, data, func, u0, l0, v0);
  cfopts.fdim = numel(f);
end
func.data = data;

fdim = cfopts.fdim;
switch type
  case { 'zero' }
    fidx                = complementary.f_dim + (1:fdim);
    complementary.f_dim = complementary.f_dim + fdim;
  case { 'inactive' 'active' 'internal' }
    fidx                = complementary.f_dim + (1:fdim);
    complementary.f_dim = complementary.f_dim + fdim;
  otherwise
    fidx                = [];
end
func.f_idx  = fidx;
complementary.f_idx = [complementary.f_idx fidx];

if ~isempty(pnames)
  assert(fdim==pnum, '%s: %s %s', mfilename, 'number of parameter', ...
    'names must match dimension of function');
end

prob = coco_add_func_after(prob, 'coco_adjoint', ...
  @complementary_add, func, args);    

complementary.funcs       = [ complementary.funcs func ];
complementary.(funcarray) = [ complementary.(funcarray) ...
  numel(complementary.funcs) ];
complementary.identifyers = [ complementary.identifyers fid ];

%% add pnames to list of parameters and update efunc.p_dim
%  check for duplicates, add names and compute parameter indices

if pnum>0
  check_names([complementary.idx2par efunc.idx2par], pnames, 'parameter');
end
complementary.idx2par = [ complementary.idx2par pnames ];

%% return efunc in prob

prob.complementary = complementary;

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
