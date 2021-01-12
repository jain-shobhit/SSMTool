function [prob, data] = ep_add_HB(prob, data, sol, varargin)
%EP_ADD_HB   Append 'ep' instance Hopf bifurcation zero problem.
%
% [PROB DATA] = EP_ADD_HB(PROB, DATA, SOL, [OPTS])
% OPTIONS = { '-no-test' }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the Hopf bifurcation zero problem for equilibrium points with
%     identifier coco_get_id(DATA.oid, 'ep.HB') and
%   - the Hopf bifurcation codim-2 test functions with identifier
%     coco_get_id(DATA.oid, 'ep.test').
%
% to the 'ep' instance with identifier coco_get_id(DATA.oid, 'ep').
%
% The corresponding function data is stored in the fields DATA.ep_hb and
% DATA.hb_tst, respectively. The instance of the 'ep' equilibrium point
% zero problem with identifier coco_get_id(DATA.oid, 'ep') must have been
% created with the option '-cache-jac' passed to EP_ADD.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'ep' toolbox instance data structure (see also
%        ODE_INIT_DATA).
% SOL  : 'ep' instance solution structure (see also EP_READ_SOLUTION). The
%        Hopf bifurcation zero problem is initialized with data in the
%        struct SOL.hb. This struct must contain the fields u0 and t0,
%        defining the initial solution guesses for each function object.
%        Typically, t0=[], which means that the continuation direction is
%        determined by the covering method.
%
% OPTS : '-no-test' (optional) disables adding of codim-2 test functions.
%
% On output:
%
% PROB : Continuation problem structure with 'ep.HB' instance added.
% DATA : Toolbox data structure with Hopf bifurcation zero problem data
%        added. EP_ADD_HB will add up to two structures to DATA: ep_hb
%        containing function data of the 'ep' Hopf bifurcation zero problem
%        and hb_tst (if test functions are added) containing additional
%        function data of the corresponding codim-2 test function. These
%        structures contain the following fields:
%
%        ep_hb.fid   : function identifier of 'ep' Hopf bifurcation zero
%            problem; fid = coco_get_id(OID, 'ep.HB');
%        ep_hb.v_idx : context-independent indices of eigenvector of J^2
%            in vector u of continuation variables; v = u(ep_hb.v_idx)
%        ep_hb.w_idx : context-independent indices of image of eigenvector
%            of J^2 in vector u of continuation variables; w =
%            u(ep_hb.w_idx)
%        ep_hb.k_idx : context-independent index of the negative of the
%            corresponding eigenvalue in vector u of continuation
%            variables; k = u(ep_hb.k_idx)
%        ep_hb.nv    : vector orthogonal to eigenvector of J^2
%        ep_hb.m     : number of rows of Jacobian of Hopf bifurcation
%            zero problem; not saved with solution data
%        ep_hb.n     : number of columns of Jacobian of Hopf bifurcation
%            zero problem; not saved with solution data
%        ep_hb.rows  : row index array used in sparse construction of
%            Jacobian of Hopf bifurcation zero problem; not saved with
%            solution data
%        ep_hb.cols  : column index array used in sparse construction of
%            Jacobian of Hopf bifurcation zero problem; not saved with
%            solution data
%        ep_hb.vals  : preallocated content for Jacobian of Hopf
%            bifurcation zero problem; not saved with solution data
%
%        hb_tst.fid  : function identifier of 'ep.HB' codim-2 test
%            functions; fid = coco_get_id(OID, 'ep.test')
%
% See also: ODE_INIT_DATA, EP_ADD, EP_READ_SOLUTION, ODE_ISOL2EP,
% ODE_EP2EP, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_add_HB.m 2902 2015-10-09 18:06:32Z hdankowicz $

opts_spec = {
    '-no-test', 'NT', false, 'toggle', {}
  };

opts = coco_parse_opts(opts_spec, varargin{:});

uidx = coco_get_func_data(prob, data.ep_var.fid, 'uidx');
gfid = coco_get_id(data.oid, 'hb_glue');
prob = coco_add_glue(prob, gfid, ...
  uidx(data.ep_var.v_idx(:,2)), uidx(data.ep_var.w_idx(:,1)));

[prob, data] = ep_HB_construct_eqn(prob, data, sol);
if ~opts.NT && data.ep.bifus
  [prob, data] = ep_HB_construct_tst(prob, data, sol);
end
end
