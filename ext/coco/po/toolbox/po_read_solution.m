function [sol, data] = po_read_solution(oid, run, varargin)
%PO_READ_SOLUTION   Read solution and toolbox data from disk.
%
% [SOL DATA] = PO_READ_SOLUTION(VARARGIN)
%
% VARARGIN = { [OID] RUN LAB }
% Read solution data from solution data file of run RUN with solution label
% LAB.
%
% On input:
%
% OID : Optional object instance identifier (string, optional).
% RUN : Run identifier (string or cell-array of strings).
% LAB : Solution label (integer).
%
% On output:
%
% SOL  : Solution structure.
% DATA : Toolbox data structure.
%
% PO_READ_SOLUTION reconstructs the solution and toolbox data structures of
% a saved periodic orbit and constructs restart information if the
% periodic orbit is a bifurcation point. More specifically, denote with
%
%   'po'     :  a branch of periodic orbits,
%   'po.SN'  :  a branch of saddle-node bifurcation periodic orbits,
%   'po.PD'  :  a branch of period-doubling bifurcation periodic orbits,
%   'po.TR'  :  a branch of torus bifurcation periodic orbits,
%
% and with 'BR(SP)' a special point detected along a branch of periodic
% orbits of type BR, for example, with 'po(TR)' a torus bifurcation point
% detected during periodic orbit continuation. DATA will always contain the
% fields of a 'po' instance. The solution structure SOL will have the
% fields
%
%   SOL.tbp : Base point times
%   SOL.xbp : Base point state vectors
%   SOL.T   : Interval duration
%   SOL.p   : Problem parameters
%
% and additional fields encoding an initial solution point as required by
% PO_ADD. Depending on the types of the solution branch and the equilibrium
% point, the return value of SOL will have the following additional fields:
%
%   'po(SN)' : For saddle-node bifurcations the structures SOL.var and
%      SOL.sn will be initialized with start data for coll_add_var and
%      po_add_SN. Specifically, for an autonomous problem, the field
%      SOL.var.v is set to a generalized eigenvector of the Jacobian of the
%      flow-T map corresponding to the eigenvalue 1, and perpendicular to
%      the vector field at the left end point. Similarly, for a
%      non-autonomous problem, the field SOL.var.v is set to an eigenvector
%      of the Jacobian of the flow-T map corresponding to the eigenvalue 1.
%
%   'po(PD)' : For period-doubling bifurcations the structures SOL.var and
%      SOL.pd will be initialized with start data for coll_add_var and
%      po_add_PD. Specifically, the field SOL.var.v is set to an
%      eigenvector of the Jacobian of the flow-T map corresponding to the
%      eigenvalue -1. In addition, the structure SOL.pd_orb is initialized
%      with data required to start continuation of a period-doubled branch
%      of periodic orbits from a period-doubling bifurcation.
%
%   'po(TR)' : For torus bifurcation points the structures SOL.var and
%      SOL.tr will be initialized with start data for coll_add_var and
%      po_add_TR. Specifically, the field SOL.var.v is set to the real and
%      imaginary parts of a complex eigenvector of the Jacobian of the
%      flow-T map corresponding to a complex eigenvalue of magnitude 1.
%
%   'po.SN' : The structures SOL.var and SOL.sn will be initialized with
%      restart data for coll_add_var and po_add_SN. Specifically, for an
%      autonomous problem, the field SOL.var.v is set to a generalized
%      eigenvector of the Jacobian of the flow-T map corresponding to the
%      eigenvalue 1, and perpendicular to the vector field at the left end
%      point. Similarly, for a non-autonomous problem, the field SOL.var.v
%      is set to an eigenvector of the Jacobian of the flow-T map
%      corresponding to the eigenvalue 1.
%
%   'po.PD' : The structures SOL.var and SOL.pd will be initialized with
%      start data for coll_add_var and po_add_PD. Specifically, the field
%      SOL.var.v is set to an eigenvector of the Jacobian of the flow-T map
%      corresponding to the eigenvalue -1. In addition, the structure
%      SOL.pd_orb is initialized with data required to start continuation
%      of a period-doubled branch of periodic orbits from a period-doubling
%      bifurcation.
%
%   'po.TR' : The structures SOL.var and SOL.tr will be initialized with
%      start data for coll_add_var and po_add_TR. Specifically, the field
%      SOL.var.v is set to the real and imaginary parts of a complex
%      eigenvector of the Jacobian of the flow-T map corresponding to a
%      complex eigenvalue of magnitude 1.
%
% See also: COCO_READ_SOLUTION, ODE_ISOL2PO, ODE_PO2PO, PO_ADD, PO_ADD_SN,
% PO_ADD_PD, PO_ADD_TR, COLL_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_read_solution.m 3148 2019-11-08 20:52:57Z hdankowicz $

if nargin<3
  [oid, run, lab] = coco_deal('', oid, run);
else
  lab = varargin{1};
end

[tbid, format, branch_type] = guess_format(oid, run, lab);
[data, chart, uidx] = coco_read_solution(tbid, run, lab, 'data', ...
  'chart', 'uidx');

sol = struct('format', format, 'branch_type', branch_type, ...
  'pt_type', chart.pt_type, 'u', chart.x, 't', chart.t);

if coco_is_chart_data(chart, 'po.test')
  sol.po_test = coco_get_chart_data(chart, 'po.test');
end

switch format
  
  case 'po.v2'
    [sol, data] = read_po_v2(sol, data, chart, uidx, tbid, branch_type, run, lab);
    
  case 'po.v1'
  %  [sol, data] = read_po_v1(sol, data, chart, uidx, tbid, branch_type);
end

end

function [tbid, format, branch_type] = guess_format(oid, run, lab)
% Guess format of solution file from contents. This function will accept
% data files from other toolboxes, if all fields required by coll are
% present.

[tb, sol_info] = coco_read_tb_info(oid, run, lab, 'tb', 'po');

if isempty(tb)
  tbid = coco_get_id(oid, 'po');
else
  tbid = coco_get_id(oid, tb);
end

if ~isempty(sol_info)
  format = sol_info.format;
  branch_type = sol_info.branch_type;
  return
else
  format = 'po.v1';
end

% Check for contents present in v1.
data = coco_read_solution(tbid, run, lab, 'data');
assert(~isempty(data), ...
  '%s: could not find solution data of toolbox instance ''%s''', ...
  mfilename, tbid);

if isfield(data, 'xbp_idx') && isfield(data, 'T_idx') ...
    && isfield(data, 'p_idx')
  branch_type = 'coll';
  return
end

error('%s: cannot restart from given solution data.', mfilename);
end

function [sol, data] = read_po_v2(sol, data, chart, uidx, tbid, ...
  branch_type, run, lab) %#ok<INUSL>

segoid = coco_get_id(tbid, 'orb');
segsol = coll_read_solution(segoid, run, lab);

sol.tbp = segsol.tbp;
sol.xbp = segsol.xbp;
sol.T   = segsol.T;
sol.p   = segsol.p;

switch branch_type
  
  case 'po'
    switch upper(sol.pt_type)
      
      case 'BP'
      case 'FP'
      case 'SN'
        tfid  = coco_get_id(tbid, 'test');
        cdata = coco_get_chart_data(chart, tfid);
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for SN continuation\n', ...
            mfilename);
        else
          sol.var.v = cdata.sn.v;
          
          sol.sn.u0 = -cdata.sn.b;
          sol.sn.t0 = [];
        end
        
      case 'PD'
        tfid  = coco_get_id(tbid, 'test');
        cdata = coco_get_chart_data(chart, tfid);
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for PD continuation\n', ...
            mfilename);
        else
          sol.var.v = cdata.pd.v;
          
          sol.pd_orb.x0 = cdata.pd.x0;
          sol.pd_orb.t0 = cdata.pd.t0;
          sol.pd_orb.p0 = cdata.pd.p0;
          
          sol.pd.u0 = [];
          sol.pd.t0 = [];
        end
        
      case 'TR'
        tfid  = coco_get_id(tbid, 'test');
        cdata = coco_get_chart_data(chart, tfid);
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for TR continuation\n', ...
            mfilename);
        else
          sol.var.v = cdata.tr.v;
          
          sol.tr.u0 = [cdata.tr.a; cdata.tr.b];
          sol.tr.t0 = [];
        end
        
    end
    
  case 'po.SN'
    fid  = coco_get_id(tbid, 'SN');
    chart = coco_read_solution(fid, run, lab, 'chart');
    sol.var.v = chart.x(data.po_sn.v_idx);
    beta = chart.x(data.po_sn.b_idx);
    sol.sn.u  = beta;
    sol.sn.t  = [];
    sol.sn.u0 = beta;
    sol.sn.t0 = [];
    
  case 'po.PD'
    fid = coco_get_id(data.po_orb.fid, 'PD');
    cdata = coco_get_chart_data(chart, fid);
    sol.pd_orb.x0 = cdata.pd.x0;
    sol.pd_orb.t0 = cdata.pd.t0;
    sol.pd_orb.p0 = cdata.pd.p0;
   
    chart = coco_read_solution(fid, run, lab, 'chart');
    sol.var.v = chart.x(data.po_pd.v_idx);
    sol.pd.u  = [];
    sol.pd.t  = [];
    sol.pd.u0 = [];
    sol.pd.t0 = [];
    
  case 'po.TR'    
    fid   = coco_get_id(tbid, 'TR');
    chart = coco_read_solution(fid, run, lab, 'chart');
    sol.var.v = chart.x(data.po_tr.v1_idx);
    a = chart.x(data.po_tr.a_idx);
    b = chart.x(data.po_tr.b_idx);
    sol.tr.u  = [a; b];
    sol.tr.t  = [];
    sol.tr.u0 = [a; b];
    sol.tr.t0 = [];

end

end
