function po = po_sol_info(sub_type)
%PO_SOL_INFO   Construct 'po' instance solution information structure.
%
% PO = PO_SOL_INFO(SUB_TYPE)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct a solution information structure appropriate for the subtype
% passed in SUB_TYPE, which may be omitted or empty to set the branch type
% to 'po'. Valid sub types are 'SN', 'PD', and 'TR', which set the branch
% type to 'po.SN', 'po.PD', and 'po.TR', respectively.  The solution
% information structure is appended to the ode toolbox information and used
% by PO_READ_SOLUTION to determine what information is stored for a
% solution point.
%
% See also: PO_READ_SOLUTION, ODE_ADD_TB_INFO

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_sol_info.m 2904 2015-10-10 02:25:50Z hdankowicz $

po.format = 'po.v2';
if nargin<1 || isempty(sub_type)
  po.branch_type = 'po';
else
  po.branch_type = [ 'po.' sub_type ];
end
end
