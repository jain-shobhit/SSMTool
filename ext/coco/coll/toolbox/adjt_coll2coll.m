function prob = adjt_coll2coll(prob, oid, varargin)
%ADJT_COLL2COLL   Append adjoint of 'coll' instance from saved solution.
%
% PROB = ADJT_COLL2COLL(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-coll-end' | '-end-coll' }
%
% Append adjoint of a 'coll' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB from the same saved solution using ODE_COLL2COLL.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of trajectory segments.
%
% See ODE_ISOL2COLL for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to restart a new continuation run.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of trajectory
%        segments. Pass non-trivial object identifiers if an instance of
%        the COLL toolbox is part of a composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% OPTS : '-switch', '-coll-end', and '-end-coll' (optional, multiple
%        options may be given). '-switch' instructs ADJT_COLL2COLL to
%        switch branches at the solution point, which must then be a branch
%        point. Either '-coll-end' or '-end-coll' marks the end of input to
%        ADJT_COLL2COLL.
%
% See also: ADJT_ISOL2COLL, ADJT_BP2COLL, COLL_READ_ADJOINT,
% COLL_ADJT_INIT_DATA, COLL_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_coll2coll.m 2898 2015-10-07 21:17:13Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
   'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
  'SOID',     '',   'str', 'soid', oid, 'read', {}
   'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-coll-end',       '',    '',    'end', {}
  '-end-coll',       '',    '',    'end', {}
    '-switch', 'switch', false, 'toggle', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

if opts.switch
  prob = adjt_BP2coll(prob, oid, args.run, args.soid, args.lab);
  return
end

tbid = coco_get_id(oid, 'coll');
data = coco_get_func_data(prob, tbid, 'data');

sol = coll_read_adjoint(args.soid, args.run, args.lab);
if ~isempty(sol.tl0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.tl0 = sol.tl;
  sol.pars_tl0 = sol.pars_tl;
  sol.T0_tl0 = sol.T0_tl;
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end

data = coll_adjt_init_data(prob, data, oid, 'coll_seg');
prob = coll_construct_adjt(prob, data, sol);
end
