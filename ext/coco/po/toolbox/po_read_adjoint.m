function [sol, data] = po_read_adjoint(oid, run, varargin)
%PO_READ_ADJOINT   Read adjoint data from disk.
%
% [SOL DATA] = PO_READ_ADJOINT(VARARGIN)
%
% VARARGIN = { [OID] RUN LAB }
% Read solution data from solution data file of run RUN with solution label
% LAB.
%
% VARARGIN = { '' '' DATA }
% Construct initial adjoint solution and data structure.
%
% On input:
%
% OID : Optional object instance identifier (string, optional).
% RUN : Run identifier (string or cell-array of strings).
% LAB : Solution label (integer).
%
% DATA : Data structure.
%
% On output:
%
% SOL  : Adjoint solution structure.
% DATA : Adjoint data structure or unchanged copy of input argument DATA.
%
% In the first calling form, PO_READ_ADJOINT reconstructs the solution and
% toolbox data structures of a saved solution and constructs
% restart information if the solution is a bifurcation point. More
% specifically, denote with
%
%   'po'     :  a branch of periodic orbits,
%
% and with 'BR(SP)' a special point detected along a branch of periodic
% orbits of type BR, for example, with 'po(PD)' a period doubling
% bifurcation point detected during periodic orbit continuation. The
% solution structure SOL will have the fields
%
%   SOL.L    :  Lagrange multipliers for boundary conditions
%   SOL.TL   :  Differentials of Lagrange multipliers for boundary
%               conditions
%
% and additional fields encoding an initial solution point as required by
% PO_CONSTRUCT_ADJT. Depending on the types of the solution branch and the
% equilibrium point, the return value of SOL will have the following
% additional fields:
%
%   'po(BP)' : For branch-points the field SOL.TL0 will be initialized to
%      a singular vector normal to SOL.TL.
%
% See also: COCO_READ_ADJOINT, ADJT_ISOL2PO, ADJT_PO2PO,
% PO_ADJT_INIT_DATA

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ep_read_solution.m 2902 2015-10-09 18:06:32Z hdankowicz $

if isempty(oid) && isempty(run)
  [sol, data] = read_sol_from(varargin{:});
  return
end

if nargin<3
  [oid, run, lab] = coco_deal('', oid, run);
else
  lab = varargin{1};
end

info = coco_read_tb_info(oid, run, lab, 'po');
format      = info.format;
branch_type = info.branch_type;

tbid = coco_get_id(oid, 'po');

switch format
  
  case 'po.v2'
    
    [data, chart1, lidx1] = coco_read_adjoint(tbid, run, lab, 'data', ...
      'chart', 'lidx');
    
    sol.l0  = chart1.x;
    sol.tl  = chart1.t;
    sol.tl0 = [];
    
    pfid    = coco_get_id(tbid, 'period');
    [chart2, lidx2] = coco_read_adjoint(pfid, run, lab, 'chart', 'lidx');
    sol.T_l0  = chart2.x;
    sol.T_tl  = chart2.t;
    sol.T_tl0 = [];
    
    if ~data.ode.autonomous
      sfid = coco_get_id(tbid, 'tinit');
      [chart3, lidx3] = coco_read_adjoint(sfid, run, lab, 'chart', 'lidx');
      sol.tinit_l0  = chart3.x;
      sol.tinit_tl  = chart3.t;
      sol.tinit_tl0 = [];
    else
      lidx3 = [];
      sol.tinit_l0 = [];
      sol.tinit_tl = [];
    end
    sol.tinit_tl0 = [];
    
    
    switch branch_type
      
      case 'po'
        switch upper(chart1.pt_type)
          
          case 'BP'
            cdata = coco_get_chart_data(chart1, 'lsol');
            if isempty(cdata)
              coco_warn([], 1, 1, ...
                '%s: could not find restart data for branch-switching\n', ...
                mfilename);
            else
              sol.tl0       = cdata.v(lidx1);
              sol.T_tl0     = cdata.v(lidx2);
              sol.tinit_tl0 = cdata.v(lidx3);
            end
        end
    end
    
  otherwise
    error('adjoint formulation not available in earlier versions of COCO');
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

sol.l0        = [];
sol.tl0       = [];
sol.T_l0      = [];
sol.T_tl0     = [];
sol.tinit_l0  = [];
sol.tinit_tl0 = [];

end
