function ep = ep_sol_info(sub_type)
%EP_SOL_INFO   Construct 'ep' instance solution information structure.
%
% EP = EP_SOL_INFO(SUB_TYPE)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct a solution information structure appropriate for the subtype
% passed in SUB_TYPE, which may be omitted or empty to set the branch type
% to 'ep'. Valid subtypes are 'VAR', 'SN', and 'HB', which set the branch
% type to 'ep.VAR', 'ep.SN', and 'ep.HB', respectively. The solution
% information structure is appended to the ode toolbox information and used
% by EP_READ_SOLUTION to determine what information is stored for a
% solution point.
%
% See also: EP_READ_SOLUTION, ODE_ADD_TB_INFO

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_sol_info.m 2897 2015-10-07 17:43:39Z hdankowicz $

ep.format = 'ep.v3';
if nargin<1 || isempty(sub_type)
  ep.branch_type = 'ep';
else
  ep.branch_type = [ 'ep.' sub_type ];
end
end
