function data = hspo_init_data(prob, src_data, oid, varargin)
%HSPO_INIT_DATA   Initialize 'hspo' instance data structure.
%
% DATA = HSPO_INIT_DATA(PROB, SRC, OID, VARARGIN)
% VARARGIN = { [NAME]... }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct 'hspo' instance data structure.
%
% Input arguments:
%
% PROB     : Continuation problem structure.
% SRC_DATA : Source data structure.
% OID      : Object identifier of 'hspo' toolbox instance.
% VARARGIN : List of toolbox field names to copy or initialize.
%
% On return, DATA contains the following fields:
%
% oid        : Object identifier, set to OID.
% hspo       : Settings of 'hspo' instance.
% bvid       : Corresponding 'bvp' instance identifier.
% fhan       : Array of function handles to RHS f of evolution equation,
%              event function h, and reset function g.
% dfdxhan    : Array of function handles to Jacobians df/dx of RHS of 
%              evolution equation, dh/dx of event function, and dgd/dx of
%              reset function.
% dfdphan    : Array of function handles to Jacobians df/dp of RHS of
%              evolution equation dh/dx of event function, and dgd/dx of
%              reset function.
% modes      : Array of mode labels.
% events     : Array of event labels.
% resets     : Array of reset labels.
% pnames     : List of parameter names.
% event_F    : Wrapper to evaluation of h.
% event_DFDX : Wrapper to evaluation of dh/dx.
% event_DFDP : Wrapper to evaluation of dh/dp.
% reset_F    : Wrapper to evaluation of g.
% reset_DFDX : Wrapper to evaluation of dg/dx.
% reset_DFDP : Wrapper to evaluation of dg/dp.
% no_save    : List of field names to be excluded by coco_save_data.
%
% and any fields with names listed in VARARGIN.
%
% The fields fhan, dfdxhan, dfdphan, modes, events, resets and pnames are
% copied from the source data structure SRC_DATA and must be present. Any
% fields with names passed in VARARGIN are either copied from SRC_DATA if a
% field with this name is present in SRC_DATA, or initialized to the empty
% structure. The function event_F and its Jacobians event_DFDX and
% event_DFDP, and similarly the function reset_F and its Jacobians
% reset_DFDX and reset_DFDP, provide wrappers to the event function and its
% Jacobians and to the jump function and its Jacobians, taking the
% presence/absence of the fields dfdxhan and dfdphan into account. The
% field no_save is initialized to the empty set and collects names of
% fields to be omitted by the slot function COCO_SAVE_DATA. The constructed
% data structure DATA is a protected instance of COCO_FUNC_DATA.
%
% See also HSPO_SETTINGS, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_init_data.m 2839 2015-03-05 17:09:01Z fschild $

data = coco_func_data('protect');
data.oid = oid;

fields = [ 'hspo' varargin ];
for i=1:numel(fields)
  field = fields{i};
  if isfield(src_data, field)
    data.(field) = src_data.(field);
  else
    data.(field) = struct();
  end
end

tbid = coco_get_id(oid, 'hspo');
data.hspo  = hspo_get_settings(prob, tbid, data.hspo);

data.fhan      = src_data.fhan;
data.dfdxhan   = src_data.dfdxhan;
data.dfdphan   = src_data.dfdphan;
data.dfdxdxhan = src_data.dfdxdxhan;
data.dfdxdphan = src_data.dfdxdphan;
data.dfdpdphan = src_data.dfdpdphan;
data.modes     = src_data.modes;
data.events    = src_data.events;
data.resets    = src_data.resets;
data.pnames    = src_data.pnames;
data.bvid      = coco_get_id(tbid, 'orb.bvp');
data.no_save   = {};

data.event_F      = @event_F;
data.event_DFDX   = @event_DFDX;
data.event_DFDP   = @event_DFDP;
data.event_DFDXDX = @event_DFDXDX;
data.event_DFDXDP = @event_DFDXDP;
data.event_DFDPDP = @event_DFDPDP;
data.reset_F      = @reset_F;
data.reset_DFDX   = @reset_DFDX;
data.reset_DFDP   = @reset_DFDP;
data.reset_DFDXDX = @reset_DFDXDX;
data.reset_DFDXDP = @reset_DFDXDP;
data.reset_DFDPDP = @reset_DFDPDP;

end

function f = event_F(data, x, p, event)
%event_F   Evaluation of event function.
f = data.fhan{2}(x, p, event);
end

function Jx = event_DFDX(data, x, p, event)
%event_DFDX   Evaluation of Jacobian of event function.

if isempty(data.dfdxhan{2})
  Jx = coco_ezDFDX('f(x,p)', data.fhan{2}, x, p, event);
else
  Jx = data.dfdxhan{2}(x, p, event);
end

end

function Jp = event_DFDP(data, x, p, event)
%event_DFDP   Evaluation of Jacobian of event function.

if isempty(data.dfdphan{2})
  Jp = coco_ezDFDP('f(x,p)', data.fhan{2}, x, p, event);
else
  Jp = data.dfdphan{2}(x, p, event);
end

end

function Jxx = event_DFDXDX(data, x, p, event)
%event_DFDXDX   Evaluation of Jacobian of event function.

if isempty(data.dfdxdxhan{2})
  Jxx = coco_ezDFDX('f(x,p)', data.dfdxhan{2}, x, p, event);
else
  Jxx = data.dfdxdxhan{2}(x, p, event);
end

end

function Jxp = event_DFDXDP(data, x, p, event)
%event_DFDXDX   Evaluation of Jacobian of event function.

if isempty(data.dfdxdphan{2})
  Jxp = coco_ezDFDP('f(x,p)', data.dfdxhan{2}, x, p, event);
else
  Jxp = data.dfdxdphan{2}(x, p, event);
end

end

function Jpp = event_DFDPDP(data, x, p, event)
%event_DFDXDX   Evaluation of Jacobian of event function.

if isempty(data.dfdpdphan{2})
  Jpp = coco_ezDFDP('f(x,p)', data.dfdphan{2}, x, p, event);
else
  Jpp = data.dfdpdphan{2}(x, p, event);
end

end

function f = reset_F(data, x, p, reset)
%reset_F   Evaluation of reset function.
f = data.fhan{3}(x, p, reset);
end

function Jx = reset_DFDX(data, x, p, reset)
%reset_DFDX   Evaluation of Jacobian of reset function.

if isempty(data.dfdxhan{3})
  Jx = coco_ezDFDX('f(x,p)', data.fhan{3}, x, p, reset);
else
  Jx = data.dfdxhan{3}(x, p, reset);
end

end

function Jp = reset_DFDP(data, x, p, reset)
%reset_DFDP   Evaluation of Jacobian of reset function.

if isempty(data.dfdphan{3})
  Jp = coco_ezDFDP('f(x,p)', data.fhan{3}, x, p, reset);
else
  Jp = data.dfdphan{3}(x, p, reset);
end

end

function Jxx = reset_DFDXDX(data, x, p, reset)
%event_DFDXDX   Evaluation of Jacobian of reset function.

if isempty(data.dfdxdxhan{3})
  Jxx = coco_ezDFDX('f(x,p)', data.dfdxhan{3}, x, p, reset);
else
  Jxx = data.dfdxdxhan{3}(x, p, reset);
end

end

function Jxp = reset_DFDXDP(data, x, p, reset)
%event_DFDXDX   Evaluation of Jacobian of event function.

if isempty(data.dfdxdphan{3})
  Jxp = coco_ezDFDP('f(x,p)', data.dfdxhan{3}, x, p, reset);
else
  Jxp = data.dfdxdphan{3}(x, p, reset);
end

end

function Jpp = reset_DFDPDP(data, x, p, reset)
%event_DFDXDX   Evaluation of Jacobian of event function.

if isempty(data.dfdpdphan{3})
  Jpp = coco_ezDFDP('f(x,p)', data.dfdphan{3}, x, p, reset);
else
  Jpp = data.dfdpdphan{3}(x, p, reset);
end

end