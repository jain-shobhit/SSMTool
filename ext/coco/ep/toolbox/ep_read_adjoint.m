function [sol, data] = ep_read_adjoint(oid, run, varargin)
%EP_READ_ADJOINT   Read adjoint data from disk.
%
% [SOL DATA] = EP_READ_ADJOINT(VARARGIN)
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
% In the first calling form, EP_READ_ADJOINT reconstructs the adjoint
% solution and data structures from a saved solution and constructs restart
% information if the solution is a branch point. More specifically, denote
% with
%
%   'ep'     :  a branch of equilibrium points
%
% and with 'BR(SP)' a special point detected along a branch of equilibrium
% points of type BR, for example, with 'ep(HB)' a Hopf bifurcation point
% detected during equilibrium point continuation. The solution structure
% SOL will have the fields
%
%   SOL.L    :  Lagrange multipliers
%   SOL.TL   :  Differentials of Lagrange multipliers
%
% and additional fields encoding an initial solution point as required by
% EP_CONSTRUCT_ADJT. Depending on the types of the solution branch and the
% equilibrium point, the return value of SOL will have the following
% additional fields:
%
%   'ep(BP)' : For branch-points the field SOL.TL0 will be initialized to
%      a singular vector normal to SOL.T.
%
% See also: COCO_READ_ADJOINT, ADJT_ISOL2EP, ADJT_EP2EP, EP_ADJT_INIT_DATA

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ep_read_adjoint.m 2902 2015-10-09 18:06:32Z hdankowicz $

if isempty(oid) && isempty(run)
  [sol, data] = read_sol_from(varargin{:});
  return
end

if nargin<3
  [oid, run, lab] = coco_deal('', oid, run);
else
  lab = varargin{1};
end

[tbid, format, branch_type] = guess_format(oid, run, lab);

switch format
  
  case 'ep.v3'
    
    [data, chart1, lidx1] = coco_read_adjoint(tbid, run, lab, 'data', ...
      'chart', 'lidx');
    
    sol.l0  = chart1.x;
    sol.tl  = chart1.t;
    sol.tl0 = [];
    
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
      
      case 'ep'
        switch upper(chart1.pt_type)
          
          case 'BP'
            cdata = coco_get_chart_data(chart1, 'lsol');
            if isempty(cdata)
              coco_warn([], 1, 1, ...
                '%s: could not find restart data for branch-switching\n', ...
                mfilename);
            else
              sol.tl0 = cdata.v(lidx1);
              sol.pars_tl0 = cdata.v(lidx2);
            end
        end
    end
    
  case {'ep.v2' 'ep.v1'}
    error('adjoint formulation not available in earlier versions of COCO');
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

sol.l0       = [];
sol.tl0      = [];
sol.pars_l0  = [];
sol.pars_tl0 = [];

end

function [tbid, format, branch_type] = guess_format(oid, run, lab)
% Guess format of solution file from contents. This function will accept
% data files from other toolboxes, if all fields required by ep are
% present.

[tb, sol_info] = coco_read_tb_info(oid, run, lab, 'tb', 'ep');

if isempty(tb)
  tbid = coco_get_id(oid, 'ep');
else
  tbid = coco_get_id(oid, tb);
end

if ~isempty(sol_info)
  format = sol_info.format;
  branch_type = sol_info.branch_type;
  return
else
  format = 'ep.v1';
end

% Check for contents present in v1.
data = coco_read_solution(tbid, run, lab, 'data');
assert(~isempty(data), ...
  '%s: could not find solution data of toolbox instance ''%s''', ...
  mfilename, tbid);

if isfield(data, 'x_idx') ...
    && isfield(data, 'p_idx')
  branch_type = 'ep';
  
  if isfield(data, 'w') ...
      && isfield(data, 'v_idx') ...
      && isfield(data, 'k_idx')
    branch_type = 'ep.HB';
  elseif isfield(data, 'v_idx')
    branch_type = 'ep.SN';
  end
  
  return
end

error('%s: cannot restart from given solution data.', mfilename);

end
