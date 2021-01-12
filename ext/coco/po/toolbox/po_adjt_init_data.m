function data = po_adjt_init_data(prob, src_data, oid, varargin) %#ok<INUSL>
%PO_ADJT_INIT_DATA   Initialize data structure for a 'po' adjoint problem.
%
% DATA = PO_ADJT_INIT_DATA(PROB, SRC, OID, VARARGIN)
% VARARGIN = { [NAME]... }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Initialize data structure for the 'po' adjoint problem. Content of the
% adjoint data structure is initially inherited from the 'po' instance
% data structure.
%
% Input arguments:
%
% PROB     : Continuation problem structure.
% SRC_DATA : Source data structure.
% OID      : Object identifier of 'po' toolbox instance.
% VARARGIN : List of toolbox field names to copy or initialize.
%
% On return, DATA contains the following fields:
%
% oid       : Object identifier, set to OID.
% cid       : 'coll' instance identifier.
% no_save   : List of field names to be excluded by coco_save_data.
%
% and any fields with names listed in VARARGIN.
%
% The field cid is copied from the source data structure SRC_DATA and must
% be present. Any fields with names passed in VARARGIN are either copied
% from SRC_DATA if a field with this name is present in SRC_DATA, or
% initialized to the empty structure. The field no_save is initialized to
% the empty set and collects names of fields to be omitted by the slot
% function COCO_SAVE_DATA. The constructed data structure DATA is a
% protected instance of COCO_FUNC_DATA.
%
% See also PO_INIT_DATA, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_init_data.m 2839 2015-03-05 17:09:01Z fschild $

data = coco_func_data('protect');
data.oid = oid;

fields = ['ode' varargin];
for i=1:numel(fields)
  field = fields{i};
  if isfield(src_data, field)
    data.(field) = src_data.(field);
  else
    data.(field) = struct();
  end
end

data.cid     = src_data.cid;
data.no_save = {};

end
