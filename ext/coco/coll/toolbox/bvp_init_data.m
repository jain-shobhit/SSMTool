function data = bvp_init_data(prob, src_data, oid, varargin) %#ok<INUSL>
%BVP_INIT_DATA   Initialize 'bvp' instance data structure.
%
% DATA = BVP_INIT_DATA(PROB, SRC, OID, VARARGIN)
% VARARGIN = { [NAME]... }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct 'bvp' instance data structure and initialize with definition of
% boundary conditions.
%
% Input arguments:
%
% PROB     : Continuation problem structure.
% SRC_DATA : Source data structure.
% OID      : Object identifier of 'bvp' toolbox instance.
% VARARGIN : List of toolbox field names to copy or initialize.
%
% On return, DATA contains the following fields:
%
% oid       : Object identifier, set to OID.
% pnames    : List of parameter names.
% fhan      : Function handle to boundary conditions.
% dfdxhan   : Function handle to Jacobian of boundary conditions.
% dfdxdxhan : Function handle to evaluation of second derivatives of
%             boundary conditions.
% bc_data   : Boundary conditions data structure.
% bc_update : Function handle to boundary conditions data structure update.
% nsegs     : Number of trajectory segments.
% no_save   : List of field names to be excluded by coco_save_data.
%
% and any fields with names listed in VARARGIN.
%
% The fields pnames, fhan, dfdxhan, bc_data, bc_update, nsegs, and basemode
% are copied from the source data structure SRC_DATA and must be present.
% Any fields with names passed in VARARGIN are either copied from SRC_DATA
% if a field with this name is present in SRC_DATA, or initialized to the
% empty structure. The field no_save is initialized to the empty set and
% collects names of fields to be omitted by the slot function
% COCO_SAVE_DATA. The constructed data structure DATA is a protected
% instance of COCO_FUNC_DATA.
%
% See also ODE_SETTINGS, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_init_data.m 2839 2015-03-05 17:09:01Z fschild $

data = coco_func_data('protect');
data.oid = oid;

fields = varargin;
for i=1:numel(fields)
  field = fields{i};
  if isfield(src_data, field)
    data.(field) = src_data.(field);
  else
    data.(field) = struct();
  end
end

data.pnames    = src_data.pnames;
data.fhan      = src_data.fhan;
data.dfdxhan   = src_data.dfdxhan;
data.dfdxdxhan = src_data.dfdxdxhan;
data.bc_data   = src_data.bc_data;
data.bc_update = src_data.bc_update;
data.nsegs     = src_data.nsegs;
data.basemode  = src_data.basemode;

data.no_save = {};

end
