function prob = coco_add_tb_info(prob, oid, varargin)
%COCO_ADD_TB_INFO   Add toolbox info to data files.
%
% PROB = COCO_ADD_TB_INFO(PROB, OID, VARARGIN)
% Add a small piece of information about a toolbox to all solution and
% bd-files of a run. Only one instance of toolbox information can added to
% any object identifyer. Toolbox information is stored as a Matlab struct.
%
% Toolbox information allows access to information about a run or solution
% file without knowing the name of the toolbox that was added under an
% object identifyer. This is useful, for example, for forwarding calls to
% specialised restart parsers within toolbox families.
%
% VARARGIN = { TB_INFO }
% Add struct TB_INFO as toolbox info.
%
% VARARGIN = INPUT_TO_STRUCT
% Add struct constructed from arguments in INPUT_TO_STRUCT as toolbox info.
%
% See also: STRUCT.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_tb_info.m 2839 2015-03-05 17:09:01Z fschild $

if nargin==3
  tb_info = varargin{1};
  assert(isstruct(tb_info) && isscalar(tb_info), ...
    '%s: toolbox info is not a struct.', mfilename);
else
  tb_info = struct();
  N = numel(varargin);
  i = 1;
  while i<=N
    tb_info.(varargin{i}) = varargin{i+1};
    i = i+2;
  end
end
fid  = coco_get_id(oid, 'tb_info');
prob = coco_add_slot(prob, fid, @coco_save_data, tb_info, 'save_reduced');
prob = coco_add_slot(prob, fid, @coco_save_data, tb_info, 'save_bd');
end
