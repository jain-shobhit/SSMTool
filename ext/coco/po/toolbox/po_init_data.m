function data = po_init_data(prob, src_data, oid, varargin)
%PO_INIT_DATA   Initialize 'po' instance data structure.
%
% DATA = PO_INIT_DATA(PROB, SRC, OID, VARARGIN)
% VARARGIN = { [NAME]... }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct 'po' instance data structure. 
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
% oid      : Object identifier, set to OID.
% po       : Settings of 'po' instance.
% cid      : Corresponding 'coll' instance identifier
% no_save  : List of field names to be excluded by coco_save_data.
%
% and any fields with names listed in VARARGIN.
%
% Any fields with names passed in VARARGIN are either copied from SRC_DATA
% if a field with this name is present in SRC_DATA, or initialized to the
% empty structure. The field no_save is initialized to the empty set and
% collects names of fields to be omitted by the slot function
% COCO_SAVE_DATA. The constructed data structure DATA is a protected
% instance of COCO_FUNC_DATA.
%
% See also PO_SETTINGS, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_init_data.m 2839 2015-03-05 17:09:01Z fschild $

data = coco_func_data('protect');
data.oid = oid;

fields = [ 'po' varargin ];
for i=1:numel(fields)
  field = fields{i};
  if isfield(src_data, field)
    data.(field) = src_data.(field);
  else
    data.(field) = struct();
  end
end

tbid = coco_get_id(oid, 'po');
data.po  = po_get_settings(prob, tbid, data.po);
data.cid = coco_get_id(tbid, 'orb.coll');
data.no_save = {};

end
