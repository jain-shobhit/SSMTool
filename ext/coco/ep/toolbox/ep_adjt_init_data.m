function data = ep_adjt_init_data(prob, src_data, oid, varargin) %#ok<INUSL>
%EP_ADJT_INIT_DATA   Initialize data structure for an 'ep' adjoint problem.
%
% DATA = EP_ADJT_INIT_DATA(PROB, SRC_DATA, OID, VARARGIN)
% VARARGIN = { [NAME]... }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Initialize data structure for the 'ep' adjoint problem. Content of the
% adjoint data structure is initially inherited from the ode family data
% structure associated with the corresponding 'ep' zero problem.
%
% Input arguments:
%
% PROB     : Continuation problem structure.
% SRC_DATA : Source data structure.
% OID      : Object identifier of ODE toolbox instance.
% VARARGIN : List of toolbox field names to copy or initialize.
%
% On return, DATA contains the following fields:
%
% oid        : Object identifier, set to OID.
% fhan       : Function handle to RHS f of evolution equation.
% dfdxhan    : Function handle to Jacobian df/dx.
% dfdphan    : Function handle to Jacobian df/dp.
% dfdxdxhan  : Function handle to Jacobian d2f/dx^2.
% dfdxdphan  : Function handle to Jacobian d2f/dxdp.
% dfdpdphan  : Function handle to Jacobian d2f/dp^2.
% xdim       : Number of state variables.
% pdim       : Number of problem parameters.
% pnames     : List of parameter names.
% ode        : Settings of ODE toolbox instance.
% ode_F      : Wrapper to vectorized evaluation of f.
% ode_DFDX   : Wrapper to vectorized evaluation of df/dx.
% ode_DFDP   : Wrapper to vectorized evaluation of df/dp.
% ode_DFDXDX : Wrapper to vectorized evaluation of d2f/dx^2.
% ode_DFDXDP : Wrapper to vectorized evaluation of d2f/dxdp.
% ode_DFDPDP : Wrapper to vectorized evaluation of d2f/dp^2.
% no_save    : List of field names to be excluded by coco_save_data.
%
% and any fields with names listed in VARARGIN.
%
% The fields fhan, dfdxhan, dfdphan, dfdxdxhan, dfdxdphan, dfdpdphan, xdim,
% pdim, and pnames are copied from the source data structure SRC_DATA and
% must be present. Any fields with names passed in VARARGIN are either
% copied from SRC_DATA, if a field with this name is present in SRC_DATA,
% or initialized to the empty structure. The fields ode_F, ode_DFDX,
% ode_DFDP, ode_DFDXDX, ode_DFDXDP, and ode_DFDPDP are initialized to
% function handles that provide wrappers for vectorized evaluation of the
% right-hand side of a given ODE and its derivatives, taking the value of
% the ODE setting 'vectorized' and the presence/absence of the fields
% dfdxhan, dfdphan, dfdxdxhan, dfdxdphan, and dfdpdphan into account. The
% field no_save is initialized to the empty set and collects names of
% fields to be omitted by the slot function COCO_SAVE_DATA. The constructed
% data structure DATA is a protected instance of COCO_FUNC_DATA.
%
% See also ODE_INIT_DATA, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ep_adjt_init_data.m 2902 2015-10-09 18:06:32Z hdankowicz $

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

data.xdim       = src_data.xdim;
data.pdim       = src_data.pdim;
data.pnames     = src_data.pnames;
data.fhan       = src_data.fhan;
data.dfdxhan    = src_data.dfdxhan;
data.dfdphan    = src_data.dfdphan;
data.dfdxdxhan  = src_data.dfdxdxhan;
data.dfdxdphan  = src_data.dfdxdphan;
data.dfdpdphan  = src_data.dfdpdphan;
data.ode_F      = src_data.ode_F;
data.ode_DFDX   = src_data.ode_DFDX;
data.ode_DFDP   = src_data.ode_DFDP;
data.ode_DFDXDX = src_data.ode_DFDXDX;
data.ode_DFDXDP = src_data.ode_DFDXDP;
data.ode_DFDPDP = src_data.ode_DFDPDP;

data.no_save = {};

end
