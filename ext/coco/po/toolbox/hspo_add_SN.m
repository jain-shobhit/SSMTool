function [prob, data] = hspo_add_SN(prob, data, sol, varargin)
%HSPO_ADD_SN   Add 'hspo' instance saddle-node bifurcation zero problem.
%
% [PROB DATA] = HSPO_ADD_SN(PROB, DATA, SOL, [OPTS])
% OPTIONS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the saddle-node bifurcation zero problem for multisegment periodic
%     orbits with identifier coco_get_id(DATA.oid, 'hspo.SN')
%
% to the 'hspo' instance with identifier coco_get_id(DATA.oid, 'hspo').
%
% The corresponding function data is stored in the field DATA.hspo_sn.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'hspo' toolbox instance data structure (see also 
%        HSPO_INIT_DATA, and ODE_ISOL2HSPO).
% SOL  : 'hspo' instance solution structure (see also HSPO_READ_SOLUTION).
%        The saddle-node bifurcation zero problem is initialized with data
%        in the struct SOL.sn. This struct must contain the fields u0 and
%        t0, defining the initial solution guesses for each function
%        object. Typically, t0=[], which means that the continuation
%        direction is determined by the covering method.
%
% OPTS : No options implemented in the current version.
%
% On output: 
%
% PROB : Continuation problem structure with 'hspo.SN' instance added.
% DATA : Toolbox data structure with saddle-node bifurcation zero problem
%        data added. HSPO_ADD_SN will add one structure to DATA: hspo_sn
%        containing function data of the 'hspo' saddle-node bifurcation
%        zero problem. This structure contains the following fields:
%
%        hspo_sn.fid   : function identifier of 'hspo' saddle-node
%            bifurcation zero problem; fid = coco_get_id(OID, 'hspo.SN');
%        hspo_sn.p_idx : context-independent indices of problem parameters
%            in vector u of continuation variables; p = u(hspo_sn.p_idx)
%        hspo_sn.x_idx : context-independent indices of collection of
%            initial end points on trajectory segments in vector u of
%            continuation variables; x = u(hspo_sn.x_idx)
%        hspo_sn.v_idx : context-independent indices of eigenvector in
%            vector u of continuation variables; v = u(hspo_sn.v_idx)
%        hspo_sn.w_idx : context-independent indices of image of
%            eigenvector in vector u of continuation variables;
%            w = u(hspo_sn.w_idx)
%        hspo_sn.nsegs : number of trajectory segments
%        hspo_sn.cids  : array of 'coll' instance identifiers
%        hspo_sn.vids  : array of 'coll.var' instance identifiers
%        hspo_sn.sidx  : array of cumulative index arrays of state-space
%            dimensions.
%        hspo_sn.cont  : permutation matrix for iterated eigenvector
%            problem.
%        hspo_sn.zrs1 and hspo_sn.zrs2 : preallocated content for Jacobian
%            of saddle-node bifurcation zero problem; not saved with
%            solution data
%
% See also: HSPO_INIT_DATA, HSPO_ADD, HSPO_ADD_PD, HSPO_ADD_TR,
% HSPO_READ_SOLUTION, ODE_ISOL2HSPO, ODE_HSPO2HSPO, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add_PD.m 2839 2015-03-05 17:09:01Z fschild $

[prob, data] = hspo_SN_construct_eqn(prob, data, sol);
end
