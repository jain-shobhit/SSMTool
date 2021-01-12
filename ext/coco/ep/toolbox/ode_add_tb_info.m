function prob = ode_add_tb_info(prob, oid, tb, sol_type, varargin)
%ODE_ADD_TB_INFO   Add ODE toolbox family information.
%
% PROB = ODE_ADD_TB_INFO(PROB, OID, TB, SOL_TYPE, VARARGIN)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Append toolbox information for the ode toolbox family to the node OID.
% The ode toolbox information is a structure with fields
%
%        fam : 'ode'
%         tb : TB
%   sol_type : SOL_TYPE
%
% and any fields passed in VARARGIN to COCO_ADD_TB_INFO.
%
% SOL_TYPE is a string that will be used by various ode toolbox family
% functions to forward calls to generic functions to the appropriate
% toolbox function.

% The VARARGIN argument must follow the argument syntax for the STRUCT
% constructor function.
%
% As an example, in the assignment
%
%    prob = ode_add_tb_info(prob, oid, 'ep', 'ep', 'ep', ep_sol_info());
%
% the toolbox identifer is 'ep', the solution type is 'ep', and additional
% solution information is stored in a field with name 'ep' and value given
% by the output of ep_sol_info. The read_solution function of the toolbox
% must be able to interpret this value, which should contain only small
% pieces of data to limit memory consumption.
%
% See also: COCO_ADD_TB_INFO, STRUCT

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_add_tb_info.m 2897 2015-10-07 17:43:39Z hdankowicz $

prob = coco_add_tb_info(prob, oid, 'fam', 'ode', 'tb', tb, ...
  'sol_type', sol_type, varargin{:});
end
