function prob = ode_BP2hspo(prob, oid, varargin)
%ODE_BP2HSPO   Switch to secondary branch of multisegment periodic orbits at branch point.
%
% PROB = ODE_BP2HSPO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-hspo-end' | '-end-hspo' | '-var' VECS }
%
% Start a continuation of multisegment periodic orbits along a secondary
% branch intersecting a previously computed branch with name RUN in a
% branch point. To start from a saved branch point, at least the name RUN
% of the continuation run and the solution label LAB must be given. The
% label LAB must be the label of a branch point.
%
% The arguments and their meaning are identical to ODE_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string).
% RUN  : Run identifier (string or cell-array of strings).
% SOID : Source object instance identifier (string, optional).
% LAB  : Solution label (integer).
%
% OPTS : '-hspo-end', '-end-hspo', and '-var' VECS (optional, multiple
%        options may be given). Either '-po-end' or '-end-po' mark the end
%        of input to ODE_BP2HSPO. The option '-var' indicates the inclusion
%        of the variational problem for the corresponding trajectory
%        segment, where the initial solution guess for the perturbations to
%        the initial condition of the orbit is given by the content of
%        VECS.
% 
% See also: ODE_HSPO2HSPO, HSPO_READ_SOLUTION, HSPO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_BP2po.m 2867 2015-07-28 01:15:22Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-hspo-end',     '', '',  'end', {}
  '-end-hspo',     '', '',  'end', {}
       '-var', 'vecs', [], 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

stbid = coco_get_id(args.soid, 'hspo');
data  = coco_read_solution(stbid, args.run, args.lab, 'data');
data  = hspo_init_data(prob, data, oid, 'hspo_orb');
if data.hspo.bifus
    prob = coco_set(prob, coco_get_id(oid, 'coll'), 'var', true);
end
tsid = coco_get_id(oid, 'hspo.orb');
ssid = coco_get_id(stbid, 'orb');
prob = ode_BP2bvp(prob, tsid, args.run, ssid, args.lab, '-var', opts.vecs);
prob = hspo_add(prob, data);
prob = ode_add_tb_info(prob, oid, 'hspo', 'hspo', 'hspo', hspo_sol_info());

end
