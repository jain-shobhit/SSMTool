function prob = ode_isol2po(prob, oid, varargin)
%ODE_ISOL2PO   Start continuation of periodic orbits from initial guess.
%
% PROB = ODE_ISOL2PO(PROB, OID, VARARGIN)
% VARARGIN = { COLL [OPTS] }
% OPTS     = { '-po-end' | '-end-po'}
%
% Start a continuation of periodic orbits of evolution equations of the
% form x'=f(t,x,p) (non-autonomous) or x'=f(x,p) (autonomous), where f is
% some non-linear function.
%
% The function f can be encoded in-line as in
%
%   func = @(t,x,p) p(1,:)-x.*(p(2,:)-x.^2);
%
% or in a separate m-file as in
%
%   function y = func(t,x,p)
%     y = p(1,:)-x.*(p(2,:)-x.^2);
%   end
%
% Note that function evaluation is vectorized by default, which can be
% disabled. For large problems it is recommended to encode explicit
% derivatives.
%
% Note that the vector field f is autonomous by default. For non-autonomous
% vector fields, use coco_set(prob, 'ode', 'autonomous', false).
%
% On input:
%
% PROB   : Continuation problem structure. This structure is either
%          automatically created and passed on by COCO when using the PO
%          toolbox implicitly, or created with COCO_PROB prior to calling
%          ODE_ISOL2PO explicitly. See below for implicit and explicit use
%          of the PO toolbox.
%
%          The continuation problem structure may contain settings defined
%          with COCO_SET, which will influence the behavior of the PO
%          toolbox and the 'coll' instance embedded in a 'po' instance.
%          Execute PO_SETTINGS on the Matlab command line to see an
%          overview of PO toolbox settings specific to single-segment
%          periodic orbits.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'po'). Pass the empty
%          string '' for a simple continuation of periodic orbits. Pass a
%          non-trivial object identifier if an instance of the PO toolbox
%          is part of a composite continuation problem.
%
% COLL   : Argument sequence for construction of a 'coll' instance.
%
% OPTS   : '-po-end' and '-end-po' (optional). Either '-po-end' or
%          '-end-po' marks the end of input to ODE_ISOL2PO. 
%
% Implicit use of PO toolbox
% --------------------------
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'po', PO_ARGS, CONT_ARGS);
%
% will call the function ODE_ISOL2PO implicitly, where the arguments
% PO_ARGS following the string 'po' will be passed to ODE_ISOL2PO.
% Any remaining arguments CONT_ARGS will be passed to the continuation
% method. When using the default 1-dimensional continuation method,
% CONT_ARGS defines the list of continuation parameters and computational
% domains.
%
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'po', PO_ARGS, CONT_ARGS);
%
% differs from the above in that the existing continuation problem
% structure prob is used in the call to ODE_ISOL2PO. Non-default values for
% PO and COLL toolbox settings may be assigned to prob before this call.
%
% Explicit use of PO toolbox
% --------------------------
% A calling sequence like
%
%   prob = coco_prob();
%   ...
%   prob = ode_isol2po(prob, '', PO_ARGS);
%   ...
%   coco(prob, RUN, [], CONT_ARGS);
%
% uses the function ODE_ISOL2PO explicitly to add a 'po' instance to the
% continuation problem structure prob. This calling sequence is equivalent
% to the implicit use of the PO toolbox, but allows for inclusion of
% user-defined monitor functions between the calls to ODE_ISOL2PO and COCO
% as demonstrated in the demos.
%
% See also: COCO_SET, PO_SETTINGS, PO_READ_SOLUTION, PO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_isol2po.m 2928 2015-10-30 14:18:59Z hdankowicz $

data = po_init_data(prob, [], oid); % initialize 'po' instance data structure
if data.po.bifus
    prob = coco_set(prob, data.cid, 'var', true);
end
str  = coco_stream(varargin{:});
tsid = coco_get_id(oid, 'po.orb');
prob = ode_isol2coll(prob, tsid, str);

opts_spec = {
	'-po-end',	'',  '',	'end',	{}
	'-end-po',  '',  '',  'end',  {}
  };
coco_parse_opts(opts_spec, str);

fdata = coco_get_func_data(prob, data.cid, 'data');
data.ode = fdata.ode;
if data.ode.autonomous
  prob = po_add(prob, data);
else
  prob = po_add(prob, data, '-no-phase');
end
prob = ode_add_tb_info(prob, oid, 'po', 'po', 'po', po_sol_info());

end
