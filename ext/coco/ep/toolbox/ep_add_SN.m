function [prob, data] = ep_add_SN(prob, data, sol, varargin)
%EP_ADD_SN   Add 'ep' instance saddle-node bifurcation zero problem.
%
% [PROB DATA] = EP_ADD_SN(PROB, DATA, SOL, [OPTS])
% OPTS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the saddle-node bifurcation zero problem for equilibrium points with
%     identifier coco_get_id(DATA.oid, 'ep.SN').
%
% to the 'ep' instance with identifier coco_get_id(DATA.oid, 'ep').
%
% The corresponding function data is stored in the field DATA.ep_sn. The
% instance of the 'ep' equilibrium point zero problem with identifier
% coco_get_id(DATA.oid, 'ep') must have been created with the option
% '-cache-jac' passed to EP_ADD.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'ep' toolbox instance data structure (see also
%        ODE_INIT_DATA).
% SOL  : 'ep' instance solution structure (see also EP_READ_SOLUTION). The
%        saddle-node bifurcation zero problem is initialized with data in
%        the struct SOL.sn. This struct must contain the fields u0 and t0,
%        defining the initial solution guesses for each function object.
%        Typically, t0=[], which means that the continuation direction is
%        determined by the covering method.
%
% OPTS : No options implemented in the current version.
%
% On output:
%
% PROB : Continuation problem structure with 'ep.SN' instance added.
% DATA : Toolbox data structure with Hopf bifurcation zero problem data
%        added. EP_ADD_SN will add one structure to DATA: ep_sn containing
%        function data of the 'ep' saddle-node bifurcation zero problem.
%        This structure contains the following fields:
%
%        ep_sn.fid   : function identifier of 'ep' saddle-node bifurcation
%            zero problem; fid = coco_get_id(OID, 'ep.SN');
%        ep_sn.v_idx : context-independent indices of eigenvector of J in
%            vector u of continuation variables; v = u(ep_sn.v_idx)
%        ep_sn.w_idx : context-independent indices of image of eigenvector
%            of J in vector u of continuation variables; w = u(ep_sn.w_idx)
%        ep_sn.m     : number of rows of Jacobian of saddle-node
%            bifurcation zero problem; not saved with solution data
%        ep_sn.n     : number of columns of Jacobian of saddle-node
%            bifurcation zero problem; not saved with solution data
%        ep_sn.rows  : row index array used in sparse construction of
%            Jacobian of saddle-node bifurcation zero problem; not saved
%            with solution data
%        ep_sn.cols  : column index array used in sparse construction of
%            Jacobian of saddle-node bifurcation zero problem; not saved
%            with solution data
%        ep_sn.vals  : preallocated content for Jacobian of saddle-node
%            bifurcation zero problem; not saved with solution data
%
% See also: ODE_INIT_DATA, EP_ADD, EP_READ_SOLUTION, ODE_ISOL2EP,
% ODE_EP2EP, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_add_SN.m 2902 2015-10-09 18:06:32Z hdankowicz $

[prob, data] = ep_SN_construct_eqn(prob, data, sol);
end
