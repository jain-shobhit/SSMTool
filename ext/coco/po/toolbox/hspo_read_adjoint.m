function [sol, data] = hspo_read_adjoint(oid, run, varargin)
%HSPO_READ_ADJOINT   Read adjoint data from disk.
%
% [SOL DATA] = HSPO_READ_ADJOINT(VARARGIN)
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
% In the first calling form, HSPO_READ_ADJOINT reconstructs the solution
% and toolbox data structures of a saved solution and constructs restart
% information if the solution is a bifurcation point. More specifically,
% denote with
%
%   'hspo'     :  a branch of periodic orbits,
%
% and with 'BR(SP)' a special point detected along a branch of periodic
% orbits of type BR, for example, with 'hspo(PD)' a period doubling
% bifurcation point detected during periodic orbit continuation. The
% solution structure SOL will have the fields
%
%   SOL.L    :  Lagrange multipliers for period
%   SOL.TL   :  Differentials of Lagrange multipliers for period
%
% and additional fields encoding an initial solution point as required by
% HSPO_CONSTRUCT_ADJT. Depending on the types of the solution branch and
% the equilibrium point, the return value of SOL will have the following
% additional fields:
%
%   'hspo(BP)' : For branch-points the field SOL.TL0 will be initialized to
%      a singular vector normal to SOL.TL.
%
% See also: COCO_READ_ADJOINT, ADJT_ISOL2HSPO, ADJT_HSPO2HSPO,
% HSPO_ADJT_INIT_DATA

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

info = coco_read_tb_info(oid, run, lab, 'hspo');
format      = info.format;
branch_type = info.branch_type;

tbid = coco_get_id(oid, 'hspo');

switch format
  
  case 'hspo.v2'
    
    data = coco_read_solution(tbid, run, lab, 'data');
    pfid = coco_get_id(tbid, 'period');
    [chart1, lidx1] = coco_read_adjoint(pfid, run, lab, 'chart', 'lidx');
    sol.period_l0  = chart1.x;
    sol.period_tl  = chart1.t;
    sol.period_tl0 = [];  
    
    switch branch_type
      
      case 'hspo'
        switch upper(chart1.pt_type)
          
          case 'BP'
            cdata = coco_get_chart_data(chart1, 'lsol');
            if isempty(cdata)
              coco_warn([], 1, 1, ...
                '%s: could not find restart data for branch-switching\n', ...
                mfilename);
            else
              sol.period_tl = cdata.v(lidx1);
            end
        end
    end
    
  otherwise
    error('adjoint formulation not available in earlier versions of COCO');
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

sol.period_l0  = [];
sol.period_tl0 = [];

end
