function [sol, data] = bvp_read_adjoint(oid, run, varargin)
%BVP_READ_ADJOINT   Read adjoint data from disk.
%
% [SOL DATA] = BVP_READ_ADJOINT(VARARGIN)
%
% VARARGIN = { [OID] RUN LAB }
% Read adjoint data from solution data file of run RUN with solution label
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
% In the first calling form, BVP_READ_ADJOINT reconstructs the solution and
% toolbox data structures of a saved solution and constructs
% restart information if the solution is a bifurcation point. More
% specifically, denote with
%
%   'segs'   :  a branch of constrained trajectory segments,
%
% and with 'BR(SP)' a special point detected along a branch of constrained
% trajectory segments of type BR. The solution structure SOL will have the
% fields
%
%   SOL.L    :  Lagrange multipliers for boundary conditions
%   SOL.TL   :  Differentials of Lagrange multipliers for boundary
%               conditions
%
% and additional fields encoding an initial solution point as required by
% BVP_CONSTRUCT_ADJT. Depending on the types of the solution branch, the
% return value of SOL will have the following additional fields:
%
%   'segs(BP)' : For branch-points the field SOL.TL0 will be initialized to
%      a singular vector normal to SOL.TL.
%
% See also: COCO_READ_ADJOINT, ADJT_ISOL2BVP, ADJT_BVP2BVP,
% BVP_ADJT_INIT_DATA

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

info = coco_read_tb_info(oid, run, lab, 'bvp');
format      = info.format;
branch_type = info.branch_type;

tbid = coco_get_id(oid, 'bvp');

switch format
  
  case 'bvp.v1'
    
    [data, chart1, lidx1] = coco_read_adjoint(tbid, run, lab, 'data', ...
      'chart', 'lidx');
    
    sol.l0  = chart1.x;
    sol.tl  = chart1.t;
    sol.tl0 = [];
    
    for i=2:data.nsegs
      shid = sprintf('shared%d', i-1);
      sfid = coco_get_id(tbid, shid);
      [chart_s, lidx_s] = coco_read_adjoint(sfid, run, lab, 'chart', 'lidx');
      sol.([shid, '_l0'])  = chart_s.x;
      sol.([shid, '_t'])   = chart_s.t;
      sol.([shid, '_tl0']) = [];
    end
    
    if ~isempty(data.pnames)
      [chart2, lidx2] = coco_read_adjoint(coco_get_id(tbid, 'pars'), run, ...
        lab, 'chart', 'lidx');
      sol.pars_l0 = chart2.x;
      sol.pars_tl = chart2.t;
    else
      lidx2 = [];
      sol.pars_l0 = [];
      sol.pars_tl = [];
    end
    sol.pars_tl0 = [];
    
    
    switch branch_type
      
      case 'segs'
        switch upper(chart1.pt_type)
          
          case 'BP'
            cdata = coco_get_chart_data(chart1, 'lsol');
            if isempty(cdata)
              coco_warn([], 1, 1, ...
                '%s: could not find restart data for branch-switching\n', ...
                mfilename);
            else
              sol.tl0      = cdata.v(lidx1);
              for i=2:data.nsegs
                shid = sprintf('shared%d', i-1);
                sol.([shid, '_tl0']) = cdata.v(lidx_s);
              end
              sol.pars_tl0 = cdata.v(lidx2);
            end
        end
    end
    
  otherwise
    error('adjoint formulation not available in earlier versions of COCO');
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

sol.l0       = [];
sol.tl0      = [];
for i=2:data.nsegs
  sfid = sprintf('shared%d', i-1);
  sol.([sfid, '_l0'])  = [];
  sol.([sfid, '_tl0']) = [];
end
sol.pars_l0  = [];
sol.pars_tl0 = [];

end
