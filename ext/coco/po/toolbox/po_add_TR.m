function [prob, data] = po_add_TR(prob, data, sol, varargin)
%PO_ADD_TR   Add 'po' instance torus bifurcation zero problem.
%
% [PROB DATA] = PO_ADD_TR(PROB, DATA, SOL, [OPTS])
% OPTIONS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the torus bifurcation zero problem for periodic orbits with
%     identifier coco_get_id(DATA.oid, 'po.TR')
%
% to the 'po' instance with identifier coco_get_id(DATA.oid, 'po').
%
% The corresponding function data is stored in the field DATA.po_tr.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'po' toolbox instance data structure (see also 
%        PO_INIT_DATA).
% SOL  : 'po' instance solution structure (see also PO_READ_SOLUTION). The
%        torus bifurcation zero problem is initialized with data in the
%        struct SOL.tr. This struct must contain the fields u0 and t0,
%        defining the initial solution guesses for each function object.
%        Typically, t0=[], which means that the continuation direction is
%        determined by the covering method.
%
% OPTS : No options implemented in the current version.
%
% On output:
%
% PROB : Continuation problem structure with 'po.TR' instance added.
% DATA : Toolbox data structure with torus bifurcation zero problem data
%        added. PO_ADD_TR will add one structure to DATA: po_tr containing
%        function data of the 'po' torus bifurcation zero problem. This
%        structure contains the following fields:
%
%        po_tr.fid   : function identifier of 'po' torus bifurcation zero
%            problem; fid = coco_get_id(OID, 'po.TR');
%        po_tr.v1_idx and po_tr.w1_idx : context-independent indices of
%            real and imaginary parts of eigenvector in vector u of
%            continuation variables; v1 = u(po_tr.v1_idx) and v2 =
%            u(po_tr.v2_idx)
%        po_tr.w1_idx and po_tr.w2_idx : context-independent indices of
%            real and imaginary parts of image of eigenvector in vector u
%            of continuation variables; w1 = u(po_tr.w1_idx) and w2 =
%            u(po_tr.w2_idx)
%        po_tr.a_idx and po_tr.b_idx : context-independent indices of
%            cosine and sine of angle of rotation in real eigenvector
%            conditions; a = u(po_tr.a_idx) and b = u(po_tr.b_idx)
%        po_tr.xbp_idx, po_pd.T_idx, po_pd.p_idx, and po_pd.var_idx :
%            context-independent indices of variational collocation problem
%            unknowns in vector u of continuation variables; not saved with
%            solution data
%        po_tr.maps and po_pd.mesh : copies of coll_seg.maps and
%            coll_seg.mesh fields of corresponding 'coll' instance; not
%            saved with solution data
%        po_tr.m     : number of rows of Jacobian of torus bifurcation 
%            zero problem; not saved with solution data
%        po_tr.n     : number of columns of Jacobian of torus bifurcation 
%            zero problem; not saved with solution data
%        po_tr.rows  : row index array used in sparse construction of
%            Jacobian of torus bifurcation zero problem; not saved with
%            solution data
%        po_tr.cols  : column index array used in sparse construction of
%            Jacobian of torus bifurcation zero problem; not saved with
%            solution data
%        po_tr.vals  : preallocated content for Jacobian of torus
%            bifurcation zero problem; not saved with solution data
%
% See also: PO_INIT_DATA, PO_ADD, PO_ADD_SN, PO_ADD_PD, PO_READ_SOLUTION, 
% ODE_ISOL2PO, ODE_PO2PO, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add_TR.m 2906 2015-10-13 19:10:23Z hdankowicz $

[prob, data] = po_TR_construct_eqn(prob, data, sol);
end
