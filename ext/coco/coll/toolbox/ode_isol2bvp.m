function prob = ode_isol2bvp(prob, oid, varargin)
%ODE_ISOL2BVP   Start continuation of collections of constrained trajectory segments from initial guess.
%
% PROB     = ODE_ISOL2BVP(PROB, OID, VARARGIN)
% VARARGIN = { ( COLL | {{COLL} ...} ) [PNAMES] BCND [OPTS] }
% BCND     = BC [DBCDX [DBXDXDX]] [BC_DATA [BC_UPDATE]]
% OPTS = { '-bvp-end' | '-end-bvp' | 'F+DF' }
%
% Start a continuation of collections of constrained trajectory segments of
% evolution equations of the form x'=f(t,x,p) (non-autonomous) or x'=f(x,p)
% (autonomous), where f is some non-linear function. with boundary
% conditions bc(data,T,x0,x1,p)=0 imposed on the trajectory end points and
% the interval durations in terms of the problem parameters. At least a
% function handle BC to a function evaluating the boundary conditions must
% be provided to start a continuation of collections of constrained
% trajectory segments.
%
% As a simple example, the sequence of commands
%
% coll_fun = @(x,p) [x(2,:); -x(2,:)-p(1,:).*x(1,:)+cos(x(3,:)); ones(1,size(x,2))];
% [t0, x0] = ode45(@(t,x) func(x,1), [0 2*pi], [0; 1; 0]);
% prob = coco_prob();
% coll_args = { coll_fun, t0, x0, 1 };
% bvp_fun  = @(~, T,x0,x1,p) [T-2*pi; x1(1:2)-x0(1:2); x1(3)-x0(3)-2*pi];
% bvp_args = { {coll_args}, 'p', bvp_fun };
% prob = ode_isol2bvp(prob, '', bvp_args{:});
%
% builds a continuation problem with dimensional deficit equal to zero,
% whose solution is a closed trajectory segment with period 2*pi. Please
% see the demos shipped with COLL for more examples of use.
%
% On input:
%
% PROB       : Continuation problem structure. This structure is either
%              automatically created and passed on by COCO when using the
%              COLL toolbox implicitly, or created with COCO_PROB prior to
%              calling COLL_ISOL2BVP explicitly. See below for implicit and
%              explicit use of the COLL toolbox.
%
%              The problem structure may contain settings defined with
%              COCO_SET, which will influence the behavior of each 'coll'
%              instance embedded in a 'bvp' instance. Execute COLL_SETTINGS
%              on the Matlab command line to see an overview of COLL
%              toolbox settings.
%
% OID        : Object instance identifier (string). The corresponding
%              toolbox instance identifier is coco_get_id(OID, 'bvp'). Pass
%              the empty string '' for a simple continuation of collections
%              of constrained trajectory segments. Pass a non-trivial
%              object identifier if an instance of the COLL toolbox is part
%              of a composite continuation problem.
%
% COLL       : Argument sequence for construction of a 'coll' instance,
%              excluding PNAMES.
%
% PNAMES     : Cell array of string labels for continuation parameters
%              (optional). These string labels will be mapped onto the
%              problem parameters as lab{1}->p(1), lab{2}->p(2), ... If
%              string labels are given, their number must match the number
%              of problem parameters. If only one problem parameter is
%              present, a string may be given instead of a cell array of
%              strings. Omitting this argument is only relevant for
%              embedded problems. In general, this argument must be given.
% BC         : Function handle to boundary conditions.
% DBCDX      : Function handle to Jacobian w.r.t. T, x0, x1, p (optional,
%              may be empty).
% DBCDXDX    : Function handle to evaluation of second derivatives w.r.t.
%              T, x0, x1, p (optional, may be empty).
% BC_DATA    : Boundary condition function data (optional, may be empty).
% BC_UPDATE  : Function handle to update function for boundary condition
%              function data (optional, may be empty).
%
% OPTS       : '-bvp-end', '-end-bvp', or 'F+DF' (optional, multiple
%               options may be given). Either '-bvp-end' or '-end-bvp'
%               marks the end of input to ODE_ISOL2BVP. The 'F+DF' option
%               indicates the simultaneous computation of boundary
%               conditions and their Jacobian by BC.
%
% Implicit use of COLL toolbox
% --------------------------
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'bvp', BVP_ARGS, CONT_ARGS);
%
% will call the function ODE_ISOL2BVP implicitly, where the arguments
% BVP_ARGS following the string 'bvp' will be passed to ODE_ISOL2BVP.
% Any remaining arguments CONT_ARGS will be passed to the continuation
% method. When using the default 1-dimensional continuation method,
% CONT_ARGS defines the list of continuation parameters and computational
% domains.
%
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'bvp', COLL_ARGS, CONT_ARGS);
%
% differs from the above in that the existing continuation problem
% structure prob is used in the call to ODE_ISOL2BVP. Non-default values
% for 'coll' instance toolbox settings may be assigned to prob before this
% call.
%
% Explicit use of COLL toolbox
% --------------------------
% A calling sequence like
%
%   prob = coco_prob();
%   ...
%   prob = ode_isolbvp(prob, '', BVP_ARGS);
%   ...
%   coco(prob, RUN, [], CONT_ARGS);
%
% uses the function ODE_ISOL2BVP explicitly to add a 'bvp' instance of the
% COLL toolbox to the problem structure. This calling sequence is
% equivalent to the implicit use of the COLL toolbox. However, it allows
% for inclusion of user-defined monitor functions between the calls to
% ODE_ISOL2BVP and COCO, as demonstrated in the demos.
%
% See also: COCO_SET, COLL_SETTINGS, BVP_READ_SOLUTION, ODE_ISOL2COLL

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_isol2bvp.m 2839 2015-03-05 17:09:01Z fschild $

tbid = coco_get_id(oid, 'bvp');
str  = coco_stream(varargin{:});
if isa(str.peek, 'function_handle')
  nsegs   = 1;
  segoid  = coco_get_id(tbid, 'seg1');
  cids{1} = coco_get_id(segoid, 'coll');
  prob    = ode_isol2coll(prob, segoid, str);
elseif isa(str.peek, 'cell')
  coll  = str.get;
  nsegs = numel(coll);
  assert(nsegs>0, '%s: insufficient number of segments', tbid);
  cids  = cell(1,nsegs);
  for i=1:nsegs
    assert(isa(coll{i}, 'cell'), '%s: incompatible argument', tbid);
    segoid  = coco_get_id(tbid, sprintf('seg%d', i));
    cids{i} = coco_get_id(segoid, 'coll');
    prob    = ode_isol2coll(prob, segoid, coll{i}{:});
  end
end

grammar   = '[PNAMES] BC [DBCDX [DBCDXDX]] [BC_DATA [BC_UPDATE]] [OPTS]';
args_spec = {
   'PNAMES', 'cell', '{str}',    'pnames',       {}, 'read', {}
       'BC',     '',     '@',      'fhan',       [], 'read', {}
    'DBCDX',     '',  '@|[]',   'dfdxhan',       [], 'read', {}
  'DBCDXDX',     '',  '@|[]', 'dfdxdxhan',       [], 'read', {}
  'BC_DATA',     '',  'stct',   'bc_data', struct(), 'read', {}
'BC_UPDATE',     '',  '@|[]', 'bc_update',       [], 'read', {}
  };
opts_spec = {
      'F+DF', 'fdf', false, 'toggle', {}
  '-bvp-end',    '',    '',    'end', {}
  '-end-bvp',    '',    '',    'end', {}
  };
[data, opts] = coco_parse(grammar, args_spec, opts_spec, str);

pnum = [];
for i=1:nsegs
  fdata = coco_get_func_data(prob, cids{i}, 'data');
  pdim  = fdata.pdim;
  assert(isempty(fdata.pnames), ...
    '%s: parameter labels must not be passed to coll', tbid);
  assert(isempty(pnum) || pnum==pdim, '%s: %s', tbid, ...
    'number of parameters must be equal for all segments');
  pnum = pdim;
end
assert(pdim==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''pnames''', tbid);

if opts.fdf
  data.basemode = 1;
elseif ~isempty(data.dfdxhan)
  data.basemode = 2;
else
  data.basemode = 3;
end
data.nsegs = nsegs;
data = bvp_init_data(prob, data, oid); % initialize 'bvp' instance data structure
data.cids  = cids;
prob = bvp_close_segs(prob, data);
prob = ode_add_tb_info(prob, oid, 'bvp', 'segs', 'bvp', bvp_sol_info());

end
