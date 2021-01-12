function prob = adjt_isol2ep(prob, oid, varargin)
%ADJT_ISOL2EP   Append adjoint of 'ep' instance from initial guess.
%
% PROB = ADJT_ISOL2EP(PROB, OID, VARARGIN)
% VARARGIN = { [OPTS] }
% OPTS = { '-ep-end' | '-end-ep' }
%
% Append adjoint of an 'ep' instance with object instance identifier OID
% that has been previously added to the continuation problem contained in
% PROB using ODE_ISOL2EP. The preceding call to ODE_ISOL2EP must include
% explicit Jacobians, while functions evaluating second derivatives are
% optional.
%
% On input:
%
% PROB : Continuation problem structure.
%
% OID  : Object instance identifier (string). The corresponding toolbox
%        instance identifier is coco_get_id(OID, 'ep'). Pass the empty
%        string '' for a simple continuation of equilibrium points. Pass a
%        non-trivial object identifier if an instance of the EP toolbox is
%        part of a composite continuation problem.
%
% OPTS : '-ep-end' and '-end-ep'(optional, multiple options may be given).
%        Either '-ep-end' or '-end-ep' mark the end of input to
%        ADJT_ISOL2EP.
%
% See also: ODE_ISOL2EP, EP_READ_ADJOINT, EP_ADJT_INIT_DATA,
% @EP_CONSTRUCT_ADJT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: adjt_isol2ep.m 2901 2015-10-09 02:47:22Z hdankowicz $

grammar   = '[OPTS]';
args_spec = {};
opts_spec = {
  '-ep-end',       '',    '',    'end', {}
  '-end-ep',       '',    '',    'end', {}
  };
coco_parse(grammar, args_spec, opts_spec, varargin{:});

tbid  = coco_get_id(oid, 'ep');
fdata = coco_get_func_data(prob, tbid, 'data');

[sol, data] = ep_read_adjoint('', '', fdata);
data = ep_adjt_init_data(prob, data, oid);
prob = ep_construct_adjt(prob, data, sol);

end
