function [prob, data] = po_add_PD(prob, data, sol, varargin)
%PO_ADD_PD   Add 'po' instance period-doubling bifurcation zero problem.
%
% [PROB DATA] = PO_ADD_PD(PROB, DATA, SOL, [OPTS])
% OPTIONS = {}
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% This function appends
%
%   - the period-doubling bifurcation zero problem for periodic orbits with
%     identifier coco_get_id(DATA.oid, 'po.PD')
%
% to the 'po' instance with identifier coco_get_id(DATA.oid, 'po').
%
% The corresponding function data is stored in the field DATA.po_pd.
%
% On input:
%
% PROB : Continuation problem structure.
% DATA : An initialized 'po' toolbox instance data structure (see also 
%        PO_INIT_DATA).
% SOL  : 'po' instance solution structure (see also PO_READ_SOLUTION). The
%        period-doubling bifurcation zero problem is initialized with data
%        in the struct SOL.pd. This struct must contain the fields u0 and
%        t0, defining the initial solution guesses for each function
%        object. Typically, t0=[], which means that the continuation
%        direction is determined by the covering method.
%
% OPTS : No options implemented in the current version.
%
% On output:
%
% PROB : Continuation problem structure with 'po.PD' instance added.
% DATA : Toolbox data structure with period-doubling bifurcation zero
%        problem data added. PO_ADD_PD will add one structure to DATA:
%        po_pd containing function data of the 'po' period-doubling
%        bifurcation zero problem. This structure contains the following
%        fields:
%
%        po_pd.fid   : function identifier of 'po' period-doubling
%            bifurcation zero problem; fid = coco_get_id(OID, 'po.PD');
%        po_pd.v_idx : context-independent indices of eigenvector in vector
%            u of continuation variables; v = u(po_pd.v_idx)
%        po_pd.w_idx : context-independent indices of image of eigenvector
%            in vector u of continuation variables; w = u(po_pd.w_idx)
%        po_pd.xbp_idx, po_pd.T_idx, po_pd.p_idx, and po_pd.var_idx :
%            context-independent indices of variational collocation problem
%            unknowns in vector u of continuation variables; not saved with
%            solution data
%        po_pd.maps and po_pd.mesh : copies of coll_seg.maps and
%            coll_seg.mesh fields of corresponding 'coll' instance; not
%            saved with solution data
%        po_pd.m     : number of rows of Jacobian of period-doubling
%            bifurcation zero problem; not saved with solution data
%        po_pd.n     : number of columns of Jacobian of period-doubling
%            bifurcation zero problem; not saved with solution data
%        po_pd.rows  : row index array used in sparse construction of
%            Jacobian of period-doubling bifurcation zero problem; not
%            saved with solution data
%        po_pd.cols  : column index array used in sparse construction of
%            Jacobian of period-doubling bifurcation zero problem; not
%            saved with solution data
%        po_pd.vals  : preallocated content for Jacobian of period-doubling
%            bifurcation zero problem; not saved with solution data
%
% See also: PO_INIT_DATA, PO_ADD, PO_ADD_SN, PO_ADD_TR, PO_READ_SOLUTION,
% ODE_ISOL2PO, ODE_PO2PO, COCO_GET_ID

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_add_PD.m 2906 2015-10-13 19:10:23Z hdankowicz $

[prob, data] = po_PD_construct_eqn(prob, data, sol);
end
