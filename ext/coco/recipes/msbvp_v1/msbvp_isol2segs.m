function prob = msbvp_isol2segs(prob, oid, varargin)
%MSBVP_ISOL2SEGS   Append 'msbvp' instance constructed from initial data.
%
% Construct multiple instances of 'coll' and append boundary conditions.
%
% PROB     = MSBVP_ISOL2SEGS(PROB, OID, VARARGIN)
% VARARGIN = { COLL... ( PNAMES | 'END-COLL' ) BCND }
% BCND     = @BC [@DBCDX] [BC_DATA [@BC_UPDATE]]
%
% PROB       - Continuation problem structure.
% OID        - Object instance identifier (string).
% COLL       - Argument to coll_isol2seg without PNAMES.
% PNAMES     - Optional string label or cell array of string labels for
%              continuation parameters tracking problem parameters.
% @BC        - Boundary conditions function handle.
% @DBCDX     - Optional function handle to Jacobian w.r.t. T, x0, x1, p.
% BC_DATA    - Optional boundary condition function data (struct).
% @BC_UPDATE - Optional function handle to function data updater.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_isol2segs.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'msbvp'); % Create toolbox instance identifier
str  = coco_stream(varargin{:});  % Convert varargin to stream of tokens for argument parsing
data.nsegs = 0;
while isa(str.peek, 'function_handle')
  data.nsegs = data.nsegs+1;
  segoid = coco_get_id(tbid, sprintf('seg%d', data.nsegs)); % Create unique object instance identifier
  prob   = coll_isol2seg(prob, segoid, str); % Use 'coll' constructor to parse one instance of COLL
end
data.pnames = {};
if strcmpi(str.peek, 'end-coll') % Check for stop token
  str.skip;
else
  data.pnames = str.get('cell');
end
data.fbchan = str.get;
data.dfbcdxhan = [];
if is_empty_or_func(str.peek)
  data.dfbcdxhan = str.get;
end
data.bc_data   = struct();
data.bc_update = [];
if isstruct(str.peek)
  data.bc_data = str.get;  
  if is_empty_or_func(str.peek)
    data.bc_update = str.get;
  end
end

msbvp_arg_check(prob, tbid, data);         % Validate input
data = msbvp_init_data(prob, tbid, data);  % Build toolbox data
prob = msbvp_close_segs(prob, tbid, data); % Append boundary conditions

end

function flag = is_empty_or_func(x)
%IS_EMPTY_OR_FUNC   Check if input is empty or contains a function handle.
flag = isempty(x) || isa(x, 'function_handle');
end
