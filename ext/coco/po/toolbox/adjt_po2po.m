function prob = adjt_po2po(prob, oid, varargin)
%ADJT_PO2PO   Append adjoint of 'po' instance from saved solution.
%
% PROB     = ADJT_PO2PO(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-po-end' | '-end-po' }
%
% Append adjoint of a 'po' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB from the same saved solution using ODE_PO2PO.
%
% PROB : Continuation problem structure.
% OID  : Target object instance identifier (string). Pass the empty string
%        '' for a simple continuation of trajectory segments.
%
% See ODE_ISOL2PO for more details on PROB and OID.
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
% OPTS : '-switch', '-po-end', and '-end-po' (optional, multiple options
%        may be given). '-switch' instructs ADJT_PO2PO to switch branches
%        at the solution point, which must then be a branch point. Either
%        '-po-end' or '-end-po' marks the end of input to ADJT_PO2PO.
%
% See also: ADJT_ISOL2PO, ADJT_BP2PO, PO_READ_ADJOINT,
% PO_ADJT_INIT_DATA, PO_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: msbvp_sol2segs.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'RUN [SOID] LAB [OPTS]';
args_spec = {
     'RUN', 'cell', '{str}',  'run',  {}, 'read', {}
    'SOID',     '',   'str', 'soid', oid, 'read', {}
     'LAB',     '',   'num',  'lab',  [], 'read', {}
  };
opts_spec = {
  '-po-end',       '',    '',    'end', {}
  '-end-po',       '',    '',    'end', {}
  '-switch', 'switch', false, 'toggle', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:}); 

if opts.switch
  prob = adjt_BP2po(prob, oid, args.run, args.soid, args.lab);
  return
end

[sol, data] = po_read_adjoint(args.soid, args.run, args.lab);
if ~isempty(sol.tl0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.tl0 = sol.tl;
  sol.tinit_tl0 = sol.tinit_tl;
  sol.T_tl0 = sol.T_tl;
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end

data = po_adjt_init_data(prob, data, oid, 'po_orb');

soid = coco_get_id(args.soid, 'po.orb');
toid = coco_get_id(oid, 'po.orb');
cid  = coco_get_id(toid, 'coll');
prob = adjt_coll2coll(prob, toid, args.run, soid, args.lab);
data.cid = cid;

prob = po_construct_adjt(prob, data, sol);

end
