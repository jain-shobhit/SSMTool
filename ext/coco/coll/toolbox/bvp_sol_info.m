function bvp = bvp_sol_info(sub_type)
%BVP_SOL_INFO   Construct COLL solution information structure.
%
% BVP = BVP_SOL_INFO(SUB_TYPE)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct a solution information structure appropriate for the subtype
% passed in SUB_TYPE, which may be omitted or empty to set the branch type
% to 'segs'. There are currently no valid sub types. The solution
% information structure is appended to the ode toolbox information and used
% by BVP_READ_SOLUTION to determine what information is stored for a
% solution point.
%
% See also: BVP_READ_SOLUTION, ODE_ADD_TB_INFO

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_sol_info.m 2845 2015-05-12 16:52:46Z hdankowicz $

bvp.format = 'bvp.v1';
if nargin<1 || isempty(sub_type)
  bvp.branch_type = 'segs';
else
  bvp.branch_type = [ 'segs.' sub_type ];
end
end
