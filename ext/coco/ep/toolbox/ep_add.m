function [prob, data] = ep_add(prob, data, sol, varargin)
%EP_ADD   Add 'ep' instance equilibrium point zero problem.
%
% [PROB DATA] = EP_ADD(PROB, DATA, SOL, [OPTS])
% OPTS = { '-cache-jac' | '-no-test' | '-no-bddat' | '-xid' NAME }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function adds an 'ep' instance consisting of
%
%   - an equilibrium point zero problem with identifier
%     coco_get_id(DATA.oid, 'ep'),
%   - the corresponding codim-1 test functions with identifier
%     coco_get_id(DATA.oid, 'ep.test'), and
%   - a slot function with identifier coco_get_id(DATA.oid, 'ep'), adding
%     the state variable vector and its norm to the bifurcation data
%     returned by the coco entry-point function or stored to disk.
%
% The corresponding function data is stored in the fields DATA.ep_eqn,
% DATA.ep_tst and DATA.ep_bddat, respectively. If the option '-cache-jac'
% is passed, the fields DATA.ep_Jx and DATA.ep_Jp are created as read-only
% and will contain the Jacobians corresponding to the most recent
% evaluation of the equilibrium point zero problem.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'ep' toolbox instance data structure (see also
%        ODE_INIT_DATA).
% SOL  : 'ep' instance solution structure (see also EP_READ_SOLUTION). The
%        fields SOL.u0 and and SOL.t0 must contain the initial solution
%        guess and an initial continuation direction, respectively.
%        Typically, SOL.t0=[], which means that the continuation direction
%        is determined by the covering method. A non-empty vector SOL.t0
%        must be given for branch-switching at branch points.
%
% OPTS : '-cache-jac', '-no-test', '-no-bddat' and '-xid' NAME (optional,
%        multiple options may be given). '-cache-jac' enables cacheing of
%        the derivatives, which is required when adding variational
%        equations. '-no-test' disables adding of codim-1 test functions
%        and '-no-bddat' disables the addition of the state variable vector
%        and its norm to the bifurcation data returned by the coco
%        entry-point function or stored to disk. The option '-xid' NAME
%        assigns coco_get_id(DATA.oid, NAME) and sprintf('||%s||',
%        coco_get_id(DATA.oid, NAME)) as the labels of the bifurcation data
%        columns containing the state variable vector and its norm,
%        respectively. NAME defaults to 'x'.
%
% On output:
%
% PROB : Continuation problem structure with an 'ep' instance added.
% DATA : Toolbox data structure with 'ep' equilibrium point zero problem
%        data added. EP_ADD will add up to three structures to DATA: ep_eqn
%        containing additional function data of the 'ep' equilibrium point
%        zero problem, ep_tst (if test functions are added) containing
%        additional function data of the corresponding codim-1 test
%        functions, and ep_bddat (if the state variable vector and its norm
%        are added to the bifurcation data) containing function data for
%        the slot function responding to the 'bddat' signal. These
%        structures contain the following fields:
%
%        ep_eqn.x_idx : context-independent indices of state variables in
%            vector u of continuation variables; x = u(ep_eqn.x_idx)
%        ep_eqn.p_idx : context-independent indices of problem parameters
%            in vector u of continuation variables; p = u(ep_eqn.p_idx)
%        ep_eqn.fid   : function identifier of 'ep' equilibrium point zero
%            problem; fid = coco_get_id(OID, 'ep')
%
%        ep_tst.la_idx1 and ep_tst.la_idx2 : eigenvalue expansion vectors
%            required by the Hopf bifurcations test function; only
%            non-empty if detection of Hopf bifurcations is enabled;
%            require n^2 double values storage space; not saved with
%            solution data
%        ep_tst.fid : function identifier of 'ep' codim-1 test functions;
%            fid = coco_get_id(OID, 'ep.test')
%
%        ep_bddat.xid : name set with option '-xid NAME'; defaults to 'x'
%
% See also: ODE_INIT_DATA, EP_READ_SOLUTION, ODE_ISOL2EP, ODE_EP2EP,
% COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_add.m 2948 2016-04-10 10:41:39Z fschild $

opts_spec = {
  '-cache-jac',  'CJ', false, 'toggle', {}
    '-no-test',  'NT', false, 'toggle', {}
   '-no-bddat', 'DAT', false, 'toggle', {}
        '-xid', 'xid',   'x',   'read', {}
  };

opts = coco_parse_opts(opts_spec, varargin{:});
tbid = coco_get_id(data.oid, 'ep');
data.ep = ep_get_settings(prob, tbid, data.ep);
[prob, data] = ep_construct_eqn(prob, data, sol, opts.CJ);
if data.ep.bifus && ~opts.NT
  [prob, data] = ep_construct_tst(prob, data);
end
if ~opts.DAT
  data.ep_bddat.xid = coco_get_id(data.oid, opts.xid);
  prob = coco_add_slot(prob, tbid, @bddat, data, 'bddat');
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data, res] = bddat(prob, data, command, varargin)
%BDDAT   Append EP bifurcation data to BD.

res = {};
switch command
  case 'init'
    xid   = data.ep_bddat.xid;
    nrmid = sprintf('||%s||_2', xid);
    maxid = sprintf('MAX(%s)', xid);
    minid = sprintf('MIN(%s)', xid);
    res   = { nrmid, maxid, minid, xid };
  case 'data'
    eqn   = data.ep_eqn;
    chart = varargin{1}; % Current chart
    uidx  = coco_get_func_data(prob, eqn.fid, 'uidx');
    u     = chart.x(uidx);
    x     = u(eqn.x_idx);
    res   = { norm(x,2), x, x, x };
end

end
