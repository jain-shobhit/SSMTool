function prob = bvp_isol2seg(prob, oid, varargin)
%BVP_ISOL2SEG   Append 'bvp' instance constructed from initial data.
%
% Construct an instance of 'coll' and append boundary conditions.
%
% Differs from bvp_v1 by providing support for updates to function data
% used to parameterize the execution of the boundary condition zero
% function.
%
% PROB     = BVP_ISOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { COLL @BC [@DBCDX] [BC_DATA [@BC_UPDATE]] }
%
% PROB       - Continuation problem structure.
% OID        - Object instance identifier (string).
% COLL       - Argument to coll_isol2seg.
% @BC        - Boundary conditions function handle.
% @DBCDX     - Optional function handle to Jacobian w.r.t. T, x0, x1, p.
% BC_DATA    - Optional boundary condition function data (struct).
% @BC_UPDATE - Optional function handle to function data updater.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: bvp_isol2seg.m 2839 2015-03-05 17:09:01Z fschild $

tbid   = coco_get_id(oid, 'bvp');  % Create toolbox instance identifier
segoid = coco_get_id(tbid, 'seg'); % Create segment object instance identifier
str    = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
prob   = coll_isol2seg(prob, segoid, str); % Construct 'coll' instance
data.fhan = str.get;
data.dfdxhan = [];
if is_empty_or_func(str.peek)
  data.dfdxhan = str.get;
end
data.bc_data   = struct();
data.bc_update = [];
if isstruct(str.peek)
  data.bc_data = str.get;
  if isa(str.peek, 'function_handle')
    data.bc_update = str.get;
  end
end

data = bvp_init_data(prob, tbid, data); % Build toolbox data
% The following function call differs from Recipes for Continuation, 1st
% edition, page 186 (redundant prob argument).
bvp_arg_check(tbid, data);              % Validate input
prob = bvp_close_seg(prob, tbid, data); % Append boundary conditions

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
