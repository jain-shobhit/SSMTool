function [prob, data] = po_add_SN(prob, data, sol, varargin)
%PO_ADD_SN   Add 'po' instance saddle-node bifurcation zero problem.
%
% [PROB DATA] = PO_ADD_SN(PROB, DATA, SOL, [OPTS])
% OPTIONS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the saddle-node bifurcation zero problem for periodic orbits with
%     identifier coco_get_id(DATA.oid, 'po.SN')
%
% to the 'po' instance with identifier coco_get_id(DATA.oid, 'po').
%
% The corresponding function data is stored in the field DATA.po_sn.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'po' toolbox instance data structure (see also 
%        PO_INIT_DATA).
% SOL  : 'po' instance solution structure (see also PO_READ_SOLUTION). The
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
% PROB : Continuation problem structure with 'po.SN' instance added.
% DATA : Toolbox data structure with saddle-node bifurcation zero problem
%        data added. PO_ADD_SN will add one structure to DATA: po_sn
%        containing function data of the 'po' saddle-node bifurcation zero
%        problem. This structure contains the following fields:
%
%        po_sn.fid   : function identifier of 'po' saddle-node bifurcation
%            zero problem; fid = coco_get_id(OID, 'po.SN');
%        po_sn.fhan, po_sn.dfdxhan, and po_sn.dfdphan : Function handles to
%            RHS f of evolution equation and its Jacobians df/dx and dfdp;
%            only for autonomous problem; not saved with solution data.
%        po_sn.ode   : Settings of ODE toolbox instance; only for
%            autonomous problem; not saved with solution data. 
%        po_sn.ode_F, po_sn.ode_DFDX and po_sn.ode_DFDP : Wrappers to
%            vectorized evaluation of f, df/dx, and df/dp; only for 
%            autonomous problem; not saved with solution data.
%        po_sn.x_idx : context-independent indices of initial end point on
%            trajectory segment in vector u of continuation variables; only
%            for autonomous problem; x = u(po_sn.x_idx)
%        po_sn.p_idx : context-independent indices of problem parameters in
%            vector u of continuation variables; only for autonomous
%            problem; p = u(po_sn.p_idx)
%        po_sn.v_idx : context-independent indices of (generalized)
%            eigenvector in vector u of continuation variables; 
%            v = u(po_sn.v_idx)
%        po_sn.w_idx : context-independent indices of image of
%            (generalized) eigenvector in vector u of continuation variables;
%            w = u(po_sn.w_idx)
%        po_sn.b_idx : context-independent indices of coefficient of
%            initial vector field in generalized eigenvector condition in
%            vector u of continuation variables; only for autonomous
%            problem; b = u(po_sn.b_idx)
%        po_sn.m     : number of rows of Jacobian of saddle-node
%            bifurcation zero problem; not saved with solution data
%        po_sn.n     : number of columns of Jacobian of saddle-node
%            bifurcation zero problem; not saved with solution data
%        po_sn.rows  : row index array used in sparse construction of
%            Jacobian of saddle-node bifurcation zero problem; not saved
%            with solution data
%        po_sn.cols  : column index array used in sparse construction of
%            Jacobian of saddle-node bifurcation zero problem; not saved
%            with solution data
%        po_sn.vals  : preallocated content for Jacobian of saddle-node
%            bifurcation zero problem; not saved with solution data
%
% See also: PO_INIT_DATA, PO_ADD, PO_ADD_PD, PO_ADD_TR, PO_READ_SOLUTION,
% ODE_ISOL2PO, ODE_PO2PO, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add_SN.m 2904 2015-10-10 02:25:50Z hdankowicz $

[prob, data] = po_SN_construct_eqn(prob, data, sol);
end
