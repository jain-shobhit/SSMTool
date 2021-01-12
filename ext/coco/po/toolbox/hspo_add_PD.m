function [prob, data] = hspo_add_PD(prob, data, sol, varargin)
%HSPO_ADD_PD   Add 'hspo' instance period-doubling bifurcation zero problem.
%
% [PROB DATA] = HSPO_ADD_PD(PROB, DATA, SOL, [OPTS])
% OPTIONS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the period-doubling bifurcation zero problem for multisegment
%     periodic orbits with identifier coco_get_id(DATA.oid, 'hspo.PD')
%
% to the 'hspo' instance with identifier coco_get_id(DATA.oid, 'hspo').
%
% The corresponding function data is stored in the field DATA.hspo_pd.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'hspo' toolbox instance data structure (see also 
%        HSPO_INIT_DATA).
% SOL  : 'hspo' instance solution structure (see also HSPO_READ_SOLUTION).
%        The period-doubling bifurcation zero problem is initialized with
%        data in the struct SOL.pd. This struct must contain the fields u0
%        and t0, defining the initial solution guesses for each function
%        object. Typically, t0=[], which means that the continuation
%        direction is determined by the covering method.
%
% OPTS : No options implemented in the current version.
%
% On output:
%
% PROB : Continuation problem structure with 'hspo.PD' instance added.
% DATA : Toolbox data structure with period-doubling bifurcation zero
%        problem data added. PO_ADD_PD will add one structure to DATA:
%        po_pd containing function data of the 'po' period-doubling
%        bifurcation zero problem. This structure contains the following
%        fields:
%
%        hspo_pd.fid   : function identifier of 'hspo' saddle-node
%            bifurcation zero problem; fid = coco_get_id(OID, 'hspo.PD');
%        hspo_pd.p_idx : context-independent indices of problem parameters
%            in vector u of continuation variables; p = u(hspo_sn.p_idx)
%        hspo_pd.x_idx : context-independent indices of collection of
%            initial end points on trajectory segments in vector u of
%            continuation variables; x = u(hspo_pd.x_idx)
%        hspo_pd.v_idx : context-independent indices of eigenvector in
%            vector u of continuation variables; v = u(hspo_pd.v_idx)
%        hspo_pd.w_idx : context-independent indices of image of
%            eigenvector in vector u of continuation variables;
%            w = u(hspo_sn.w_idx)
%        hspo_pd.nsegs : number of trajectory segments
%        hspo_pd.cids  : array of 'coll' instance identifiers
%        hspo_pd.vids  : array of 'coll.var' instance identifiers
%        hspo_pd.sidx  : array of cumulative index arrays of state-space
%            dimensions.
%        hspo_pd.cont  : permutation matrix for iterated eigenvector
%            problem.
%        hspo_pd.zrs1 and hspo_pd.zrs2 : preallocated content for Jacobian
%            of saddle-node bifurcation zero problem; not saved with
%            solution data
%        hspo_pd.segs  : array of information extracted from coll_seg.maps,
%            coll_seg.mesh, and coll_var fields of individual 'coll'
%            instances; not saved with solution data
%
% See also: HSPO_INIT_DATA, HSPO_ADD, HSPO_ADD_SN, HSPO_ADD_TR,
% HSPO_READ_SOLUTION, ODE_ISOL2HSPO, ODE_HSPO2HSPO, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add_PD.m 2839 2015-03-05 17:09:01Z fschild $

[prob, data] = hspo_PD_construct_eqn(prob, data, sol);
end
