function prob = ode_isol2ep(prob, oid, varargin)
%ODE_ISOL2EP   Start continuation of equilibrium points from initial guess.
%
% PROB = ODE_ISOL2EP(PROB, OID, VARARGIN)
% VARARGIN = { F [DFDX [DFDP [DFDXDX [DFDXDP [DFDPDP]]]]] X0 [PNAMES] P0 [OPTS] }
% OPTS = { '-ep-end' | '-end-ep' | '-var' VECS }
%
% Start a continuation of dynamic equilibrium points of evolution equations
% of the form x'=f(x,p), or of static equilibrium points of an equation of
% the form f(x,p) = 0, where f is some non-linear function. At least a
% function handle F to a function evaluating f, an initial guess X0 for the
% state vector x, and initial values P0 for the problem parameters p must
% be provided to start a continuation of equilibrium points.
%
% The function f can be encoded in-line as in
%
%   func = @(x,p) p(1,:)-x.*(p(2,:)-x.^2);
%
% or in a separate m-file as in
%
%   function y = func(x,p)
%     y = p(1,:)-x.*(p(2,:)-x.^2);
%   end
%
% Note that function evaluation is vectorized by default. For
% non-vectorized encodings use coco_set(prob, 'ode', 'vectorized', false).
% For large problems it is recommended to encode explicit derivatives.
%
% As a simple example, the sequence of commands
%
%   func = @(x,p) p(1,:)-x.*(p(2,:)-x.^2);
%   bd = coco('1', 'ode', 'isol', 'ep', ...
%     func, 0, {'ka' 'la'}, [0; 0.5], ...
%     {'ka' 'la'}, [-0.5 0.5]);
%
% will start a continuation with run name RUN='1' of equilibrium points of
% x'=ka-x*(la-x^2) under variation of ka within the computational domain
% -0.5<=ka<=0.5, and return the computed bifurcation data as the cell array
% bd. Please see the demos shipped with EP for more examples of use. The
% meaning of the arguments EP_ARGS = { func, 0, {'ka' 'la'}, [0; 0.5] } and
% CONT_ARGS = { {'ka' 'la'}, [-0.5 0.5] } is explained below.
%
% On input:
%
% PROB   : Continuation problem structure. This structure is either
%          automatically created and passed on by COCO when using the EP
%          toolbox implicitly, or created with COCO_PROB prior to calling
%          ODE_ISOL2EP implicitly or explicitly. See below for implicit and
%          explicit use of the EP toolbox.
%
%          The continuation problem structure may contain settings defined
%          with COCO_SET, which will influence the behavior of the EP
%          toolbox. Execute EP_SETTINGS on the Matlab command line to see
%          an overview of EP toolbox settings.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'ep'). Pass the empty
%          string '' for a simple continuation of equilibrium points. Pass
%          a non-trivial object identifier if an instance of the EP toolbox
%          is part of a composite continuation problem.
%
% F      : Function handle to vector field.
% DFDX   : Function handle to Jacobian of vector field w.r.t. state
%          variables (optional, may be empty).
% DFDP   : Function handle to Jacobian of vector field w.r.t. problem
%          parameters (optional, may be empty).
% DFDXDX : Function handle to d2f/dx^2 (optional, may be empty).
% DFDXDP : Function handle to d2f/dxdp (optional, may be empty).
% DFDPDP : Function handle to d2f/dp^2 (optional, may be empty).
%
% X0     : Initial solution guess for state variables.
% PNAMES : Cell array of string labels for continuation parameters
%          (optional). These string labels will be mapped onto the problem
%          parameters as lab{1}->p(1), lab{2}->p(2), ... If string labels
%          are given, their number must match the number of problem
%          parameters. If only one problem parameter is present, a string
%          may be given instead of a cell array of strings. Omitting this
%          argument is only relevant for embedded problems. In general,
%          this argument must be given.
% P0     : Initial solution guess for problem parameters. Pass an empty
%          array to indicate absence of problem parameters.
%
% OPTS   : '-ep-end', '-end-ep', and '-var' VECS  (optional, multiple
%          options may be given). Either '-ep-end' or '-ep-coll' mark the
%          end of input to ODE_ISOL2EP. The option '-var' indicates the
%          inclusion of the variational problem J*v=w, where the initial
%          solution guess for v is given by the content of VECS.
%
% Implicit use of EP toolbox
% --------------------------
% A call to COCO of the form
%
%   coco(RUN, 'ode', 'isol', 'ep', EP_ARGS, CONT_ARGS);
%
% will call the function ODE_ISOL2EP implicitly with an empty continuation
% problem structure and empty string object instance identifier. The
% arguments EP_ARGS following the string 'ep' will be passed as VARARGIN to
% ODE_ISOL2EP. Any remaining arguments CONT_ARGS will be passed to the
% continuation method. When using the default 1-dimensional continuation
% method, CONT_ARGS defines the list of continuation parameters and
% computational domains.
%
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'ep', EP_ARGS, CONT_ARGS);
%
% differs from the above in that the existing continuation problem
% structure prob is used in the call to ODE_ISOL2EP. Non-default values for
% EP toolbox settings may be assigned to prob before this call.
%
% Explicit use of EP toolbox
% --------------------------
% A calling sequence like
%
%   prob = coco_prob();
%   ...
%   prob = ode_isol2ep(prob, '', EP_ARGS);
%   ...
%   coco(prob, RUN, [], CONT_ARGS);
%
% uses the function ODE_ISOL2EP explicitly to add an instance of the EP
% toolbox to the continuation problem structure prob. This calling sequence
% is equivalent to the implicit use of the EP toolbox, but allows for
% inclusion of user-defined monitor functions between the calls to
% ODE_ISOL2EP and COCO, as demonstrated in the demos.
%
% See also: COCO_SET, EP_SETTINGS, EP_READ_SOLUTION, EP_ADD, EP_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_isol2ep.m 2995 2017-01-21 04:33:56Z hdankowicz $

grammar   = 'F [DFDX [DFDP [DFDXDX [DFDXDP [DFDPDP]]]]] X0 [PNAMES] P0 [OPTS]';
args_spec = {
       'F',     '',     '@',      'fhan', [], 'read', {}
    'DFDX',     '',  '@|[]',   'dfdxhan', [], 'read', {}
    'DFDP',     '',  '@|[]',   'dfdphan', [], 'read', {}
  'DFDXDX',     '',  '@|[]',   'dfdthan', [], 'read', {}
  'DFDXDP',     '',  '@|[]', 'dfdxdxhan', [], 'read', {}
  'DFDPDP',     '',  '@|[]', 'dfdxdphan', [], 'read', {}
      'X0',     '', '[num]',        'x0', [], 'read', {}
  'PNAMES', 'cell', '{str}',    'pnames', {}, 'read', {}
      'P0',     '', '[num]',        'p0', [], 'read', {}
  };
opts_spec = {
  '-ep-end',     '', '',  'end', {}
  '-end-ep',     '', '',  'end', {}
     '-var', 'vecs', [], 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

[sol, data] = ep_read_solution('', '', args);
 % ode_init_data assumes that the first derivative with respect to time is
 % defined if second derivatives are given. The content of args_spec and
 % the following code makes the appropriate reassignment to accommodate the
 % restricted syntax for the 'ep' toolbox.
data.dfdpdphan = [];
data.dfdtdxhan = [];
data.dfdtdphan = [];
data.dfdtdthan = [];
data = ode_init_data(prob, data, oid, 'ep');
[data.dfdxdxhan, data.dfdxdphan, data.dfdpdphan] = ...
  deal(data.dfdthan, data.dfdxdxhan, data.dfdxdphan);

if ~isempty(opts.vecs)
  assert(isnumeric(opts.vecs) && data.xdim == size(opts.vecs,1), ...
    '%s: incompatible specification of vectors of perturbations', ...
    mfilename);
  [prob, data] = ep_add(prob, data, sol, '-cache-jac');
  prob = ep_add_var(prob, data, opts.vecs);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info('VAR'));
else
  prob = ep_add(prob, data, sol);
  prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info());
end

end
