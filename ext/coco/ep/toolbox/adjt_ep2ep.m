function prob = adjt_ep2ep(prob, oid, varargin)
%ADJT_EP2EP   Append adjoint of 'ep' instance from saved solution.
%
% PROB = ODE_EP2EP(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB [OPTS] }
% OPTS = { '-switch' | '-ep-end' | '-end-ep' }
%
% Append adjoint of an 'ep' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB from the same saved solution using ODE_EP2EP. 
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
% OPTS : '-switch', '-ep-end', and '-end-ep' (optional, multiple options
%        may be given). '-switch' instructs ADJT_EP2EP to switch branches
%        at the solution point, which must then be a branch point. Either
%        '-ep-end' or '-end-ep' mark the end of input to ADJT_EP2EP.
%
% See also: ADJT_ISOL2EP, ADJT_BP2EP, EP_READ_ADJOINT, EP_ADJT_INIT_DATA,
% EP_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: adjt_ep2ep.m 2901 2015-10-09 02:47:22Z hdankowicz $

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
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

if opts.switch
  prob = adjt_BP2ep(prob, oid, args.run, args.soid, args.lab);
  return
end

[sol, data] = ep_read_adjoint(args.soid, args.run, args.lab);
if ~isempty(sol.tl0)
  % We restart at branch point but stay on original solution manifold =>
  % reset start direction to tangent vector of primary branch.
  sol.tl0 = sol.tl;
  sol.pars_tl0 = sol.pars_tl;
  if ~coco_exist('NullItMX', 'class_prop', prob, 'cont')
    % Compute improved tangent to new branch on same solution manifold to
    % allow change of parameters on restart at a branch-point.
    prob = coco_set(prob, 'cont', 'NullItMX', 1);
  end
end
data = ep_adjt_init_data(prob, data, oid);
prob = ep_construct_adjt(prob, data, sol);

end
