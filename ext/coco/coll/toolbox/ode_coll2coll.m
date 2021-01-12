function prob = ode_coll2coll(prob, oid, varargin)
%ODE_COLL2COLL   Start continuation of trajectory segments from saved solution.
%
% PROB = ODE_COLL2COLL(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-coll-end' | '-end-coll' | '-var' VECS }
%
% Restart a continuation of trajectory segments from a trajectory segment
% that was obtained and saved to disk in a previous continuation. To
% restart from a saved trajectory segment, at least the name RUN of the
% continuation run and the solution label LAB must be given.
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
% OPTS : '-switch', '-coll-end', '-end-coll', and '-var' VECS (optional,
%        multiple options may be given). '-switch' instructs ODE_COLL2COLL
%        to switch branches at the solution point, which must then be a
%        branch point. Either '-coll-end' or '-end-coll' marks the end of
%        input to ODE_COLL2COLL. The option '-var' indicates the inclusion
%        of the variational problem, where the initial solution guess for
%        the perturbations to the trajectory initial conditions is given by
%        the content of VECS.
%
% See also: ODE_ISOL2COLL, COLL_READ_SOLUTION, COLL_ADD, COLL_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_coll2coll.m 3032 2017-09-18 02:53:12Z hdankowicz $

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
       '-var',   'vecs',    [],   'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

if opts.switch
  prob = ode_BP2coll(prob, oid, args.run, args.soid, args.lab, ...
    '-var', opts.vecs);
  return
end

[sol, data] = coll_read_solution(args.soid, args.run, args.lab);

if ~isempty(sol.t0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.t0     = sol.t;
  sol.xbp_t0 = sol.xbp_t;
  sol.T_t0   = sol.T_t;
  sol.p_t0   = sol.p_t;
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end

data = ode_init_data(prob, data, oid, 'coll');

if ~isempty(opts.vecs)
  assert(isnumeric(opts.vecs) && data.xdim == size(opts.vecs,1), ...
    '%s: incompatible specification of vectors of perturbations', ...
    mfilename);
  [prob, data] = coll_add(prob, data, sol, '-no-var', '-cache-jac');
  prob = coll_add_var(prob, data, opts.vecs);
  prob = ode_add_tb_info(prob, oid, 'coll', 'seg', 'coll', ...
    coll_sol_info('VAR'));
elseif isfield(sol, 'var')
  [prob, data] = coll_add(prob, data, sol, '-no-var', '-cache-jac');
  prob = coll_add_var(prob, data, sol.var.v);
  prob = ode_add_tb_info(prob, oid, 'coll', 'seg', 'coll', ...
    coll_sol_info('VAR'));
else
  prob = coll_add(prob, data, sol);
  prob = ode_add_tb_info(prob, oid, 'coll', 'seg', 'coll', coll_sol_info());
end
end
