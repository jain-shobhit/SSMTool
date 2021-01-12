function hspo = hspo_sol_info(sub_type)
%HSPO_SOL_INFO   Construct 'hspo' instance solution information structure.
%
% HSPO = HSPO_SOL_INFO(SUB_TYPE)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct a solution information structure appropriate for the subtype
% passed in SUB_TYPE, which may be omitted or empty to set the branch type
% to 'hspo'. Valid sub types are 'SN', 'PD', and 'TR', which set the branch
% type to 'hspo.SN', 'hspo.PD', and 'hspo.TR', respectively.  The solution
% information structure is appended to the ode toolbox information and
% used by HSPO_READ_SOLUTION to determine what information is stored for a
% solution point.
%
% See also: HSPO_READ_SOLUTION, ODE_ADD_TB_INFO

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_sol_info.m 2839 2015-03-05 17:09:01Z fschild $

hspo.format = 'hspo.v2';
if nargin<1 || isempty(sub_type)
  hspo.branch_type = 'hspo';
else
  hspo.branch_type = [ 'hspo.' sub_type ];
end
end
