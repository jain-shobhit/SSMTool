function coll = coll_sol_info(sub_type)
%COLL_SOL_INFO   Construct COLL solution information structure.
%
% COLL = COLL_SOL_INFO(SUB_TYPE)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct a solution information structure appropriate for the subtype
% passed in SUB_TYPE, which may be omitted or empty to set the branch type
% to 'seg'. The only valid subtype is 'VAR', which sets the branch type to
% 'seg.VAR'. The solution information structure is appended to the ode
% toolbox information and used by COLL_READ_SOLUTION to determine what
% information is stored for a solution point.
%
% See also: COLL_READ_SOLUTION, ODE_ADD_TB_INFO

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_sol_info.m 2898 2015-10-07 21:17:13Z hdankowicz $

coll.format = 'coll.v1';
if nargin<1 || isempty(sub_type)
  coll.branch_type = 'seg';
else
  coll.branch_type = [ 'seg.' sub_type ];
end
end
