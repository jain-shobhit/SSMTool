function [prob, data] = hspo_add(prob, data, varargin)
%HSPO_ADD   Add 'hspo' instance monitor functions.
%
% [PROB DATA] = HSPO_ADD(PROB, DATA, [OPTS])
% OPTS = { '-no-test' }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function adds an instance of
%
%   - the 'hspo' codim-1 test functions with identifier
%     coco_get_id(DATA.oid, 'hspo.test') and
%   - the orbital period monitor function with identifier
%     coco_get_id(DATA.oid, 'hspo.period').
%
% The corresponding function data is stored in the field DATA.hspo_tst.
%
% On input:
%
% PROB : Continuation problem structure. 
% DATA : An initialized 'hspo' toolbox instance data structure (see also
%        HSPO_INIT_DATA). 
% OPTS : '-no-test' (optional). '-no-test' disables adding of codim-1 test
%        functions.
%
% On output:
%
% PROB : Continuation problem structure with an 'hspo' instance added.
% DATA : Toolbox data structure with 'hspo' test function data added.
%        HSPO_ADD adds one structure to DATA: hspo_tst (if test functions
%        are added) containing additional function data of the 'hspo'
%        codim-1 test functions. This structure contains the following
%        fields:
%
%        coll_seg.x_idx : indices of state vector in vector u of continuation
%            variables; x=u(ep_eqn.x_idx)
%        ep_eqn.p_idx : indices of parameters in vector u of continuation
%            variables; p=u(ep_eqn.p_idx)
%        ep_eqn.fid   : function identifier of EP zero problem, identical
%            to the EP toolbox identifier = coco_get_id(OID, 'ep')
%
%        ep_tst.la_idx1 and ep_tst.la_idx2 : eigenvalue expansion vectors
%            required by the Hopf test function; only non-empty if Hopf
%            detection is enabled; require n^2 double values storage space
%        ep_tst.fid : function identifier of EP test functions
%            = coco_get_id(ep_eqn.fid, 'test')
%
%        ep_bddat.xid : name set with option '-xid NAME'; defaults to 'x';
%            the bddat slot is added with slot function identifier equal to
%            the toolbox identifier ep_eqn.fid
%
% See also: PO_INIT_DATA, PO_READ_SOLUTION, ODE_ISOL2PO, ODE_PO2PO,
% COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add.m 2849 2015-05-17 20:32:46Z hdankowicz $

opts_spec = {
     '-no-test',  'NT', false, 'toggle', {}
  };
opts = coco_parse_opts(opts_spec, varargin{:});

if ~opts.NT && data.hspo.bifus
  [prob, data] = hspo_construct_tst(prob, data);
end

% Add period monitor function
[fdata, uidx]= coco_get_func_data(prob, data.bvid, 'data', 'uidx');
fid  = coco_get_id(data.oid, 'hspo');
prob = coco_add_func(prob, fid, @empty, data, 'regular', {}, ...
  'uidx', uidx, 'fdim', 0);
pfid = coco_get_id(fid, 'period');
prob = coco_add_functionals(prob, pfid, ones(1,fdata.nsegs), 0, ...
  uidx(fdata.bvp_bc.T_idx), pfid, 'active');
prob = coco_add_slot(prob, fid, @coco_save_data, data, 'save_full');

end

function [data, y] = empty(prob, data, u) %#ok<INUSD,INUSL>
y = [];
end