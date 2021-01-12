function [sol, data] = hspo_read_solution(oid, run, varargin)
%HSPO_READ_SOLUTION   Read solution and toolbox data from disk.
%
% [SOL DATA] = HSPO_READ_SOLUTION(VARARGIN)
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
% HSPO_READ_SOLUTION reconstructs the solution and toolbox data structures
% of a saved multisegment periodic orbit and constructs restart information
% if the periodic orbit is a bifurcation point. More specifically, denote
% with
%
%   'hspo'     :  a branch of periodic orbits,
%   'hspo.SN'  :  a branch of saddle-node bifurcation periodic orbits,
%   'hspo.PD'  :  a branch of period-doubling bifurcation periodic orbits,
%   'hspo.TR'  :  a branch of torus bifurcation periodic orbits,
%
% and with 'BR(SP)' a special point detected along a branch of periodic
% orbits of type BR, for example, with 'po(TR)' a torus bifurcation point
% detected during periodic orbit continuation. The solution structure SOL
% will have the fields
%
%   SOL.tbp :  Base point times (segment-by-segment array)
%   SOL.xbp :  Base point state vectors (segment-by-segment array)
%   SOL.T   :  Interval durations (segment-by-segment array)
%   SOL.p   :  Problem parameters
%
% and additional fields encoding an initial solution point as required by
% HSPO_ADD. Depending on the types of the solution branch and the
% equilibrium point, the return value of SOL will have the following
% additional fields:
%
%   'hspo(SN)' : For saddle-node bifurcations the structure SOL.var will be
%      initialized with start data for coll_add_var and hspo_add_SN.
%      Specifically, the field SOL.var.v is set to a segment-by-segment
%      array whose first element is an eigenvector of the Jacobian of the
%      flow-T map corresponding to the eigenvalue 1.
%
%   'hspo(PD)' : For period-doubling bifurcations the structures SOL.var
%      and SOL.pd will be initialized with start data for coll_add_var and
%      po_add_PD. Specifically, the field SOL.var.v is set to a
%      segment-by-segment array whose first element is an eigenvector of
%      the Jacobian of the flow-T map corresponding to the eigenvalue -1.
%      In addition, the structure SOL.pd_orb is initialized with data
%      required to start continuation of a period-doubled branch of
%      multisegment periodic orbits from a period-doubling bifurcation.
%
%   'hspo(TR)' : For torus bifurcation points the structures SOL.var and
%      SOL.tr will be initialized with start data for coll_add_var and
%      po_add_TR. Specifically, the field SOL.var.v is set to a
%      segment-by-segment array whose first element contains the real and
%      imaginary parts of a complex eigenvector of the Jacobian of the
%      flow-T map corresponding to a complex eigenvalue of magnitude 1.
%
%   'hspo.SN' : The structures SOL.var and SOL.sn will be initialized with
%      restart data for coll_add_var and hspo_add_SN. Specifically, the
%      field SOL.var.v is set to a segment-by-segment array whose first
%      element is an eigenvector of the Jacobian of the flow-T map
%      corresponding to the eigenvalue -1.
%
%   'hspo.PD' : The structures SOL.var and SOL.pd will be initialized with
%      start data for coll_add_var and hspo_add_PD. Specifically, the field
%      SOL.var.v is set to a segment-by-segment array whose first element
%      is an eigenvector of the Jacobian of the flow-T map corresponding to
%      the eigenvalue -1. In addition, the structure SOL.pd_orb is
%      initialized with data required to start continuation of a
%      period-doubled branch of periodic orbits from a period-doubling
%      bifurcation.
%
%   'hspo.TR' : The structures SOL.var and SOL.tr will be initialized with
%      start data for coll_add_var and hspo_add_TR. Specifically, the field
%      SOL.var.v is set to a segment-by-segment array whose first element
%      contains the real and imaginary parts of a complex eigenvector of
%      the Jacobian of the flow-T map corresponding to a complex eigenvalue
%      of magnitude 1.
%
% See also: COCO_READ_SOLUTION, ODE_ISOL2HSPO, ODE_HSPO2HSPO, HSPO_ADD,
% HSPO_ADD_SN, HSPO_ADD_PD, HSPO_ADD_NS, COLL_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_read_solution.m 2863 2015-07-26 22:19:05Z hdankowicz $

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

switch format
  
  case 'hspo.v2'
    [sol, data] = read_hspo_v2(sol, data, chart, uidx, tbid, branch_type, run, lab);
end

end

function [tbid, format, branch_type] = guess_format(oid, run, lab)
% Guess format of solution file from contents. This function will accept
% data files from other toolboxes, if all fields required by coll are
% present.

try
  [tb, sol_info] = coco_read_tb_info(oid, run, lab, 'tb', 'hspo');
catch
  [tb, sol_info] = coco_read_tb_info(oid, run, lab, 'tb', 'po');
end

if isempty(tb)
  tbid = coco_get_id(oid, 'hspo');
else
  tbid = coco_get_id(oid, tb);
end

if ~isempty(sol_info)
  format = sol_info.format;
  branch_type = sol_info.branch_type;
  return
else
  format = 'hspo.v1';
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

function [sol, data] = read_hspo_v2(sol, data, chart, uidx, tbid, ...
  branch_type, run, lab) %#ok<INUSL>

bvpoid = coco_get_id(tbid, 'orb');
bvpsol = bvp_read_solution(bvpoid, run, lab);

for i=1:data.hspo_orb.nsegs
  sol.tbp{i} = bvpsol{i}.tbp;
  sol.xbp{i} = bvpsol{i}.xbp;
  sol.T{i}   = bvpsol{i}.T;
end
sol.p        = bvpsol{i}.p;

switch branch_type
  
  case 'hspo'
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
          
          sol.sn.u0    = [];
          sol.sn.t0    = [];
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
          sol.pd_orb.modes  = cdata.pd.modes;
          sol.pd_orb.events = cdata.pd.events;
          sol.pd_orb.resets = cdata.pd.resets;
          
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
    
  case 'hspo.SN'
    fid = coco_get_id(tbid, 'SN');
    chart = coco_read_solution(fid, run, lab, 'chart');
    v = chart.x(data.hspo_sn.v_idx);
    off = 0;
    for i=1:data.hspo_orb.nsegs
      sol.var.v{i} = v(off + (1:data.hspo_orb.dim(i)));
      off = data.hspo_orb.dim(i);
    end
    
    sol.sn.u  = [];
    sol.sn.t  = [];
    sol.sn.u0 = [];
    sol.sn.t0 = [];
    
  case 'hspo.PD'
    fid = coco_get_id(tbid, 'PD');
    cdata = coco_get_chart_data(chart, fid);
    sol.pd_orb.x0 = cdata.pd.x0;
    sol.pd_orb.t0 = cdata.pd.t0;
    sol.pd_orb.p0 = cdata.pd.p0;
    sol.pd_orb.modes  = cdata.pd.modes;
    sol.pd_orb.events = cdata.pd.events;
    sol.pd_orb.resets = cdata.pd.resets;
    
    chart = coco_read_solution(fid, run, lab, 'chart');
    v = chart.x(data.hspo_pd.v_idx);
    off = 0;
    for i=1:data.hspo_orb.nsegs
      sol.var.v{i} = v(off + (1:data.hspo_orb.dim(i)));
      off = data.hspo_orb.dim(i);
    end
    
    sol.pd.u  = [];
    sol.pd.t  = [];
    sol.pd.u0 = [];
    sol.pd.t0 = [];
    
  case 'hspo.TR'
    fid   = coco_get_id(tbid, 'TR');
    chart = coco_read_solution(fid, run, lab, 'chart');
    v = chart.x(data.hspo_tr.v_idx);
    v = reshape(v, [data.hspo_orb.cdim 2]);
    off = 0;
    for i=1:data.hspo_orb.nsegs
      sol.var.v{i} = v(off + (1:data.hspo_orb.dim(i)),:);
      off = data.hspo_orb.dim(i);
    end
    a = chart.x(data.hspo_tr.a_idx);
    b = chart.x(data.hspo_tr.b_idx);
    
    sol.tr.u  = [a; b];
    sol.tr.t  = [];
    sol.tr.u0 = [a; b];
    sol.tr.t0 = [];
    
end

end
