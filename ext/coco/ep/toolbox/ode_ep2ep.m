function prob = ode_ep2ep(prob, oid, varargin)
%ODE_EP2EP   Start continuation of equilibrium points from saved solution.
%
% PROB = ODE_EP2EP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-ep-end' | '-end-ep' | '-var' VECS | '-no-var' }
%
% Restart a continuation of equilibrium points from an equilibrium point
% that was obtained and saved to disk in a previous continuation. To
% restart from a saved equilibrium point, at least the name RUN of the
% continuation run and the solution label LAB must be given.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of equilibrium points.
%
% See ODE_ISOL2EP for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to restart a new continuation run.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of equilibrium
%        points. Pass non-trivial object identifiers if an instance of the
%        EP toolbox is part of a composite continuation problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        equilibrium point during the continuation run RUN.
%
% OPTS : '-switch', '-ep-end', '-end-ep', '-var' VECS (optional,
%        multiple options may be given), and '-no-var'. '-switch' instructs
%        ODE_EP2EP to switch branches at the solution point, which must
%        then be a branch point. Either '-ep-end' or '-end-ep' mark the end
%        of input to ODE_EP2EP. The option '-var' indicates the inclusion
%        of the variational problem J*v=w, where the initial solution guess
%        for v is given by the content of VECS. Alternatively, the option
%        '-no-var' indicates the exclusion of the variational problem.
%
% See also: ODE_ISOL2EP, EP_READ_SOLUTION, EP_ADD, EP_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_ep2ep.m 2951 2017-01-10 14:35:52Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-ep-end',       '',    '',    'end', {}
  '-end-ep',       '',    '',    'end', {}
  '-switch', 'switch', false, 'toggle', {}
     '-var',   'vecs',    [],   'read', {}
  '-no-var',  'novar', false, 'toggle', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

if opts.switch
  prob = ode_BP2ep(prob, oid, args.run, args.soid, args.lab, ...
    '-var', opts.vecs);
  return
end

[sol, data] = ep_read_solution(args.soid, args.run, args.lab);

if ~isempty(sol.t0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.t0 = sol.t;
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end

data = ode_init_data(prob, data, oid, 'ep');

if ~isempty(opts.vecs)
  assert(isnumeric(opts.vecs) && data.xdim == size(opts.vecs,1), ...
    '%s: incompatible specification of vectors of perturbations', ...
    mfilename);
  [prob, data] = ep_add(prob, data, sol, '-cache-jac');
  prob = ep_add_var(prob, data, opts.vecs);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('VAR'));
elseif isfield(sol, 'var') && ~opts.novar
  [prob, data] = ep_add(prob, data, sol, '-cache-jac');
  prob = ep_add_var(prob, data, sol.var.v);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('VAR'));
else
  prob = ep_add(prob, data, sol);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info());
end

end
