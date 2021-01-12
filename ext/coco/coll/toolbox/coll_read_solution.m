function [sol, data] = coll_read_solution(oid, run, varargin)
%COLL_READ_SOLUTION   Read solution and toolbox data from disk.
%
% [SOL DATA] = COLL_READ_SOLUTION(VARARGIN)
%
% VARARGIN = { [OID] RUN LAB } 
% Read solution data from solution data file of run RUN with solution label
% LAB.
%
% VARARGIN = { '' '' DATA }
% Construct solution data structure from initial solution guess for
% time discretization, sampled time history of state vector, and parameter
% values from content stored in fields with names 't0', 'x0', and 'p0' of
% the structure DATA.
%
% On input:
%
% OID  : Optional object instance identifier (string, optional).
% RUN  : Run identifier (string or cell-array of strings).
% LAB  : Solution label (integer).
%
% DATA : Data structure with fields t0, x0, and p0.
%
% On output:
%
% SOL  : Solution structure. 
% DATA : Toolbox data structure or copy of input argument DATA with added
%        fields xdim (state-space dimension) and pdim (parameter-space
%        dimension).
%
% In the first calling form, COLL_READ_SOLUTION reconstructs the solution
% and toolbox data structures of a saved trajectory segment. More
% specifically, denote with
%
%   'seg'     :  a branch of trajectory segments and
%   'seg.VAR' :  a branch of trajectory segments with variational problem.
%
% DATA will be a corresponding read-only instance of COCO_FUNC_DATA and
% always contain the fields of the ODE toolbox family.
%
% In both cases, the solution structure SOL will have the fields
%
% SOL.tbp : time discretization
% SOL.xbp : sampled time history
% SOL.T   : interval duration
% SOL.p   : problem parameters
%
% and additional fields encoding an initial solution point as required by
% COLL_ADD. Depending on the types of the solution branch, the return value
% of SOL will have the following additional fields:
%
%   'seg.var' : The structure SOL.var will be initialized with restart data
%      for coll_add_var.
%
% See also: COCO_READ_SOLUTION, COLL_ISOL2SEG, COLL_SEG2SEG, COLL_ADD,
% COLL_ADD_VAR

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_read_solution.m 3032 2017-09-18 02:53:12Z hdankowicz $

if isempty(oid) && isempty(run)
  [sol, data] = read_sol_from(varargin{:});
  return
end

if nargin<3
  [oid, run, lab] = coco_deal('', oid, run);
else
  lab = varargin{1};
end

try
  info = coco_read_tb_info(oid, run, lab, 'coll');
catch
  info = coco_read_tb_info(oid, run, lab, 'seg');
end
format      = info.format;
branch_type = info.branch_type;

tbid = coco_get_id(oid, 'coll');
[data, chart, uidx] = coco_read_solution(tbid, run, lab, 'data', ...
  'chart', 'uidx');

sol = struct('format', format, 'branch_type', branch_type, ...
  'pt_type', chart.pt_type, 'u', chart.x, 't', chart.t);

switch format
  
  case 'coll.v1'
    [sol, data] = read_coll_v1(sol, data, chart, uidx, tbid, ...
      branch_type, run, lab);
    
end

end

function [sol, data] = read_sol_from(data)
% Construct solution structure from data.

t0 = data.t0;
x0 = data.x0;
p0 = data.p0;

data.xdim = size(x0,2);
data.pdim = numel(p0);

assert(ndims(t0)==2 && min(size(t0))==1, ...
  '%s: input for ''t0'' is not a vector', mfilename); %#ok<*ISMAT>
assert(ndims(x0)==2, ...
  '%s: input for ''x0'' is not an array of vectors', mfilename);
assert(size(x0,1)==numel(t0), ...
  '%s: dimensions of ''t0'' and ''x0'' do not match', mfilename);
assert(numel(p0)==numel(data.pnames) || isempty(data.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', ...
  mfilename);

sol.T0  = t0(1);
sol.T   = t0(end)-t0(1);
sol.tbp = t0(:);
sol.xbp = x0;
sol.p   = p0(:);

sol.t0  = []; % Continuation direction determined by the atlas algorithm

end

function [sol, data] = read_coll_v1(sol, data, chart, uidx, tbid, branch_type, run, lab) 

seg  = data.coll_seg;
maps = seg.maps;
mesh = seg.mesh;

sol.T0  = sol.u(maps.T0_idx);
sol.T   = sol.u(maps.T_idx);
sol.tbp = sol.T0+sol.T*mesh.tbp(maps.tbp_idx);
xbp     = reshape(sol.u(maps.xbp_idx), maps.xbp_shp)';
sol.xbp = xbp(maps.tbp_idx,:);
sol.p   = sol.u(maps.p_idx);

sol.t0  = [];

switch branch_type
  case 'seg'
    switch upper(sol.pt_type)
      
      case 'BP'
        cdata = coco_get_chart_data(chart, 'lsol');
        if isempty(cdata)
          coco_warn([], 1, 1, ...
            '%s: could not find restart data for branch-switching\n', ...
            mfilename);
        else
          xbp_t      = reshape(sol.t(maps.xbp_idx), maps.xbp_shp)';
          sol.xbp_t  = xbp_t(maps.tbp_idx,:);
          sol.T0_t = sol.t(maps.T0_idx);
          sol.T_t = sol.t(maps.T_idx);
          sol.p_t    = sol.t(maps.p_idx);
          sol.t0     = cdata.v(uidx);
          xbp_t0     = reshape(sol.t0(maps.xbp_idx), maps.xbp_shp)';
          sol.xbp_t0 = xbp_t0(maps.tbp_idx,:);
          sol.T0_t0 = sol.t0(maps.T0_idx);
          sol.T_t0 = sol.t0(maps.T_idx);
          sol.p_t0   = sol.t0(maps.p_idx);
        end
        
    end
    
  case 'seg.var'
    fid = coco_get_id(tbid, 'var');
    chart = coco_read_solution(fid, run, lab, 'chart');
    sol.var.v = chart.x(data.ep_var.v0_idx);
end
end
