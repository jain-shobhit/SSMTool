function prob = ode_BP2coll(prob, oid, varargin)
%ODE_BP2COLL   Switch to secondary branch of trajectory segments at branch point.
%
% PROB = ODE_BP2COLL(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-coll-end' | '-end-coll' | '-var' VECS }
%
% Start a continuation of trajectory segments along a secondary branch
% intersecting a previously computed branch with name RUN in a branch
% point. To start from a saved branch point, at least the name RUN of the
% continuation run and the solution label LAB must be given. The label LAB
% must be the label of a branch point.
%
% The arguments and their meaning are identical to ODE_COLL2COLL.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-coll-end', '-end-coll',  and '-var' VECS (optional, multiple
%        options may be given). Either '-coll-end' or '-end-coll' mark the
%        end of input to ODE_BP2COLL. The option '-var' indicates the
%        inclusion of the variational problem, where the initial solution
%        guess for the perturbations to the trajectory initial conditions
%        is given by the content of VECS.
% 
% See also: ODE_COLL2COLL, COLL_READ_SOLUTION, COLL_ADD, COLL_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_BP2coll.m 3032 2017-09-18 02:53:12Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-coll-end',     '', '',  'end', {}
  '-end-coll',     '', '',  'end', {}
       '-var', 'vecs', [], 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
  % Compute improved tangent to secondary branch to avoid detection of
  % spurious FP that can happen at transcritical bifurcation point.
  prob = coco_set(prob, 'cont', 'NullItMX', 1);
end

[sol, data] = coll_read_solution(args.soid, args.run, args.lab);
data = ode_init_data(prob, data, oid, 'coll');

if ~isempty(opts.vecs)
  [prob, data] = coll_add(prob, data, sol, '-no-var', '-cache-jac');
  prob = coll_add_var(prob, data, args.vecs);
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
