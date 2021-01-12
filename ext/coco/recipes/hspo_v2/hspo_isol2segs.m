function prob = hspo_isol2segs(prob, oid, varargin)
%HSPO_ISOL2SEGS   Append 'hspo' instance constructed from initial data.
%
% Construct arguments for msbvp_isol2segs with boundary conditions
% corresponding to the termination of segments on event surfaces and the
% connectivity between segments imposed by the reset function.
%
% Differs from hspo_v1 by support for toolbox settings and bifurcation
% detection.
%
% PROB     = HSPO_ISOL2SEGS(PROB, OID, VARARGIN)
% VARARGIN = { FUNCS SIG {T0...} {X0...} [PNAMES] P0 }
% FUNCS    = FUN [ DX [ DP ] ]
% FUN      = { @F @E @R }
% DX       = { (@DFDX | '[]') (@DEDX | '[]') (@DRDX | '[]') }
% DP       = { (@DFDP | '[]') (@DEDP | '[]') (@DRDP | '[]') }
% SIG      = { MODE... } { EVENT... } { RESET... }
%
% PROB   - Continuation problem structure.
% OID    - Object instance identifier (string).
% @F     - Function handle to vector field.
% @E     - Function handle to event function.
% @R     - Function handle to reset function.
% @DFDX  - Optional function handle to Jacobian of vector field w.r.t.
%          problem variables.
% @DEDX  - Optional function handle to Jacobian of event function w.r.t.
%          problem variables.
% @DRDX  - Optional function handle to Jacobian of reset function w.r.t.
%          problem variables.
% @DFDP  - Optional function handle to Jacobian of vector field w.r.t.
%          problem parameters.
% @DEDP  - Optional function handle to Jacobian of event function w.r.t.
%          problem parameters.
% @DRDP  - Optional function handle to Jacobian of reset function w.r.t.
%          problem parameters.
% MODE   - Mode identifier (string).
% EVENT  - Event identifier (string).
% RESET  - Reset identifier (string).
% T0     - Array of temporal mesh points.
% X0     - Array of state vectors at mesh points.
% PNAMES - String label or cell array of string labels for
%          continuation parameters tracking problem parameters.
% P0     - Initial solution guess for parameter values.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_isol2segs.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'hspo'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
fhan = str.get('cell');
dfdxhan = cell(1,3);
dfdphan = cell(1,3);
if is_empty_or_func(str.peek('cell'))
  dfdxhan = str.get('cell');
  if is_empty_or_func(str.peek('cell'))
    dfdphan = str.get('cell');
  end
end
modes  = str.get('cell');
events = str.get('cell');
resets = str.get('cell');
t0     = str.get('cell');
x0     = str.get('cell');
pnames = {};
if iscellstr(str.peek('cell'))
  pnames = str.get('cell');
end
p0 = str.get;

hspo_arg_check(tbid, fhan, dfdxhan, dfdphan, ...
  modes, events, resets, t0, x0, p0, pnames); % Validate input
coll = {};
for i=1:numel(modes)
  coll = [coll, {...
    coll_func(fhan{1},    modes{i}), ...
    coll_func(dfdxhan{1}, modes{i}), ...
    coll_func(dfdphan{1}, modes{i}), ...
    t0{i}, x0{i}, p0}]; % Construct sequence of 'coll' arguments for msbvp_isol2segs
end
hspo_bc_data = hspo_get_settings(prob, tbid); % Get toolbox settings
hspo_bc_data = hspo_init_data(hspo_bc_data, fhan, dfdxhan, dfdphan, ...
  modes, events, resets, x0, p0); % Build boundary conditions data
prob = msbvp_isol2segs(prob, oid, coll{:}, pnames, ...
  @hspo_bc, @hspo_bc_DFDX, hspo_bc_data); % Build 'msbvp' instance
if hspo_bc_data.hspo.bifus
  prob = hspo_add_bifus(prob, oid, tbid, hspo_bc_data); % Append monitor functions
end

end

function fhan = coll_func(fhan, mode)
%COLL_FUNC   Create 'coll'-compatible function handle.
%
% Wrapper to accommodate that 'coll' compatible vector fields depend only
% on x and p. Different vector fields are identified by the mode argument
% and handled within the functional encoding using a switch statement.

if isa(fhan, 'function_handle')
  fhan = @(x,p) fhan(x,p,mode);
end
end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if elements of input are empty or contain a function handle.
flag = all(cellfun('isempty', x) | cellfun('isclass', x, 'function_handle'));
end
