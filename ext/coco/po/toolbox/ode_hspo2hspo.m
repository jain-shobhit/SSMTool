function prob = ode_hspo2hspo(prob, oid, varargin)
%ODE_HSPO2HSPO   Start continuation of multisegment periodic orbits from saved solution.
%
% PROB = ODE_HSPO2HSPO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-hspo-end' | '-end-hspo' | '-var' VECS }
%
% Restart a continuation of multisegment periodic orbits from a
% multisegment periodic orbit that was obtained and saved to disk in a
% previous continuation. To restart from a saved periodic orbit, at least
% the name RUN of the continuation run and the solution label LAB must be
% given.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of multisegment periodic orbits.
%
% See ODE_ISOL2HSPO for more details on PROB and OID.
%
% RUN  : Run identifier (string or cell-array of strings). Name of the run
%        from which to restart a new continuation run.
% SOID : Source object instance identifier (string, optional). If the
%        argument SOID is omitted, OID is used. Pass the empty string ''
%        for OID and omit SOID for a simple continuation of multisegment
%        periodic orbits. Pass non-trivial object identifiers if an
%        instance of the PO toolbox is part of a composite continuation
%        problem.
% LAB  : Solution label (integer). The integer label assigned by COCO to an
%        trajectory segment during the continuation run RUN.
%
% OPTS : '-switch', '-hspo-end', '-end-hspo', and '-var' VECS (optional,
%        multiple options may be given). '-switch' instructs ODE_HSPO2HSPO
%        to switch branches at the solution point, which must then be a
%        branch point. Either '-hspo-end' or '-end-hspo' mark the end of
%        input to ODE_HSPO2HSPO. The option '-var' indicates the inclusion
%        of a variational problem for each trajectory segment, where the
%        initial solution guess for the perturbations to the initial
%        conditions of the i-th trajectory segment is given by the content
%        of the i-th element of the cell array VECS.
%
% See also: ODE_ISOL2HSPO, HSPO_READ_SOLUTION, HSPO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_po2po.m 2849 2015-05-17 20:32:46Z hdankowicz $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-hspo-end',       '',    '',    'end', {}
  '-end-hspo',       '',    '',    'end', {}
    '-switch', 'switch', false, 'toggle', {}
       '-var',   'vecs',    [],   'read', {}
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
if opts.switch
  prob = ode_BP2bvp(prob, tsid, args.run, ssid, args.lab, ...
    '-var', opts.vecs);
else
  prob = ode_bvp2bvp(prob, tsid, args.run, ssid, args.lab, ...
    '-var', opts.vecs);
end
prob = hspo_add(prob, data);
prob = ode_add_tb_info(prob, oid, 'hspo', 'hspo', 'hspo', hspo_sol_info());

end
