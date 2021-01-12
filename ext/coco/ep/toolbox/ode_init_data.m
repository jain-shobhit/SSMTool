function data = ode_init_data(prob, src_data, oid, varargin)
%ODE_INIT_DATA   Initialize ODE toolbox family data structure.
%
% DATA = ODE_INIT_DATA(PROB, SRC_DATA, OID, VARARGIN)
% VARARGIN = { [NAME]... }
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% Construct ODE toolbox family data structure and initialize data structure
% with definition of an ODE with right-hand side f. The constructed data
% structure is suitable for use as function data for toolboxes implementing
% computations related to ODEs. For autonomous ODEs the right-hand side f
% must be defined in a function f that expects two arguments x (state
% variables) and p (problem parameters) and returns x':=dx/dt as in
%
%   x' = f(x,p) .
%
% For non-autonomous ODEs the right-hand side f must be defined in a
% function f that expects three arguments t (time), x (state variables) and
% p (problem parameters) and returns x' as in
%
%   x' = f(t,x,p) .
%
% In addition, for non-autonomous ODEs the ODE toolbox setting autonomous
% must be set to off/false.
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
% dfdthan    : Function handle to Jacobian df/dt.
% dfdxdxhan  : Function handle to Jacobian d2f/dx2.
% dfdxdphan  : Function handle to Jacobian d2f/dxdp.
% dfdpdphan  : Function handle to Jacobian d2f/dp2.
% dfdtdxhan  : Function handle to Jacobian d2f/dtdx.
% dfdtdphan  : Function handle to Jacobian d2f/dtdp.
% dfdtdthan  : Function handle to Jacobian d2f/dt2.
% xdim       : Number of state variables.
% pdim       : Number of problem parameters.
% pnames     : List of parameter names.
% ode        : Settings of ODE toolbox instance.
% ode_F      : Wrapper to vectorized evaluation of f.
% ode_DFDX   : Wrapper to vectorized evaluation of df/dx.
% ode_DFDP   : Wrapper to vectorized evaluation of df/dp.
% ode_DFDT   : Wrapper to vectorized evaluation of df/dt.
% ode_DFDXDX : Wrapper to vectorized evaluation of d2f/dx2.
% ode_DFDXDP : Wrapper to vectorized evaluation of d2f/dxdt.
% ode_DFDPDP : Wrapper to vectorized evaluation of d2f/dp2.
% ode_DFDTDX : Wrapper to vectorized evaluation of d2f/dtdx.
% ode_DFDTDP : Wrapper to vectorized evaluation of d2f/dtdp.
% ode_DFDTDT : Wrapper to vectorized evaluation of d2f/dt2.
% no_save    : List of field names to be excluded by coco_save_data.
%
% and any fields with names listed in VARARGIN.
%
% The fields fhan, dfdxhan, dfdphan, xdim, pdim, and pnames are copied from
% the source data structure SRC_DATA and must be present. For
% non-autonomous vector fields, the field dfdthan is copied from the source
% data structure SRC_DATA and must be present. Any fields with names passed
% in VARARGIN are either copied from SRC_DATA, if a field with this name is
% present in SRC_DATA, or initialized to the empty structure. The fields
% ode_F, ode_DFDX, ode_DFDP, and ode_DFDT are initialized to function
% handles that provide wrappers for vectorized evaluation of the right-hand
% side of a given ODE and its Jacobians, taking the value of the ODE
% setting 'vectorized' and the presence/absence of the fields dfdxhan,
% dfdphan, and dfdthan into account. The field no_save is initialized to
% the empty set and collects names of fields to be omitted by the slot
% function COCO_SAVE_DATA. The constructed data structure DATA is a
% protected instance of COCO_FUNC_DATA.
%
% See also ODE_SETTINGS, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_init_data.m 2995 2017-01-21 04:33:56Z hdankowicz $

data = coco_func_data('protect');
data.oid = oid;

fields = [ 'ode' varargin ];
for i=1:numel(fields)
  field = fields{i};
  if isfield(src_data, field)
    data.(field) = src_data.(field);
  else
    data.(field) = struct();
  end
end

tbid = coco_get_id(oid, 'ode');
data.ode = ode_get_settings(prob, tbid, data.ode);

data.fhan      = src_data.fhan;
data.dfdxhan   = src_data.dfdxhan;
data.dfdphan   = src_data.dfdphan;
data.dfdthan   = src_data.dfdthan;
data.dfdxdxhan = src_data.dfdxdxhan;
data.dfdxdphan = src_data.dfdxdphan;
data.dfdpdphan = src_data.dfdpdphan;
data.dfdtdxhan = src_data.dfdtdxhan;
data.dfdtdphan = src_data.dfdtdphan;
data.dfdtdthan = src_data.dfdtdthan;
data.pnames    = src_data.pnames;
data.xdim      = src_data.xdim;
data.pdim      = src_data.pdim;
data.no_save   = {};

if data.ode.autonomous
  data.ode_F      = @aut_F;
  data.ode_DFDX   = @aut_DFDX;
  data.ode_DFDP   = @aut_DFDP;
  data.ode_DFDT   = @aut_DFDT;
  data.ode_DFDXDX = @aut_DFDXDX;
  data.ode_DFDXDP = @aut_DFDXDP;
  data.ode_DFDPDP = @aut_DFDPDP;
  data.ode_DFDTDX = @aut_DFDTDX;
  data.ode_DFDTDP = @aut_DFDTDP;
  data.ode_DFDTDT = @aut_DFDTDT;
else
  data.ode_F      = @het_F;
  data.ode_DFDX   = @het_DFDX;
  data.ode_DFDP   = @het_DFDP;
  data.ode_DFDT   = @het_DFDT;
  data.ode_DFDXDX = @het_DFDXDX;
  data.ode_DFDXDP = @het_DFDXDP;
  data.ode_DFDPDP = @het_DFDPDP;
  data.ode_DFDTDX = @het_DFDTDX;
  data.ode_DFDTDP = @het_DFDTDP;
  data.ode_DFDTDT = @het_DFDTDT;
end

end

function f = aut_F(data, ~, x, p)
%F   Vectorized evaluation of autonomous F.

if data.ode.vectorized
  f = data.fhan(x, p);
else
  [m, n] = size(x);
  f     = zeros(m,n);
  fhan  = data.fhan;
  for i=1:n
    f(:,i) = fhan(x(:,i), p(:,i));
  end
end
end

function Jx = aut_DFDX(data, ~, x, p)
%DFDX   Vectorized evaluation of autonomous dF/dx.

if data.ode.vectorized
  if isempty(data.dfdxhan)
    Jx = coco_ezDFDX('f(x,p)v', data.fhan, x, p);
  else
    Jx = data.dfdxhan(x, p);
  end
else
  if isempty(data.dfdxhan)
    Jx = coco_ezDFDX('f(x,p)', data.fhan, x, p);
  else
    [m, n] = size(x);
    Jx = zeros(m,m,n);
    dfdxhan = data.dfdxhan;
    for i=1:n
      Jx(:,:,i) = dfdxhan(x(:,i), p(:,i));
    end
  end
end
end

function Jp = aut_DFDP(data, ~, x, p)
%DFDP   Vectorized evaluation of autonomous dF/dp.

if data.ode.vectorized
  if isempty(data.dfdphan)
    Jp = coco_ezDFDP('f(x,p)v', data.fhan, x, p);
  else
    Jp = data.dfdphan(x, p);
  end
else
  if isempty(data.dfdphan)
    Jp = coco_ezDFDP('f(x,p)', data.fhan, x, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jp = zeros(m,o,n);
    dfdphan = data.dfdphan;
    for i=1:n
      Jp(:,:,i) = dfdphan(x(:,i), p(:,i));
    end
  end
end
end

function Jt = aut_DFDT(~, ~, x, ~)
%DFDT   Vectorized evaluation of non-autonomous dF/dt.
Jt = zeros(size(x));
end

function Jxx = aut_DFDXDX(data, ~, x, p)
%DFDXDX   Vectorized evaluation of autonomous d2F/dx2.

if data.ode.vectorized
  if isempty(data.dfdxdxhan)
    Jxx = coco_ezDFDX('f(x,p)v', data.dfdxhan, x, p);
  else
    Jxx = data.dfdxdxhan(x, p);
  end
else
  if isempty(data.dfdxdxhan)
    Jxx = coco_ezDFDX('f(x,p)', data.dfdxhan, x, p);
  else
    [m, n] = size(x);
    Jxx = zeros(m,m,m,n);
    dfdxdxhan = data.dfdxdxhan;
    for i=1:n
      Jxx(:,:,:,i) = dfdxdxhan(x(:,i), p(:,i));
    end
  end
end
end

function Jxp = aut_DFDXDP(data, ~, x, p)
%DFDX   Vectorized evaluation of autonomous d2F/dxdp.

if data.ode.vectorized
  if isempty(data.dfdxdphan)
    Jxp = coco_ezDFDP('f(x,p)v', data.dfdxhan, x, p);
  else
    Jxp = data.dfdxdphan(x, p);
  end
else
  if isempty(data.dfdxdphan)
    Jxp = coco_ezDFDP('f(x,p)', data.dfdxhan, x, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jxp = zeros(m,m,o,n);
    dfdxdphan = data.dfdxdphan;
    for i=1:n
      Jxp(:,:,:,i) = dfdxdphan(x(:,i), p(:,i));
    end
  end
end
end

function Jpp = aut_DFDPDP(data, ~, x, p)
%DFDX   Vectorized evaluation of autonomous d2F/dp2.

if data.ode.vectorized
  if isempty(data.dfdpdphan)
    Jpp = coco_ezDFDP('f(x,p)v', data.dfdphan, x, p);
  else
    Jpp = data.dfdpdphan(x, p);
  end
else
  if isempty(data.dfdpdphan)
    Jpp = coco_ezDFDP('f(x,p)', data.dfdphan, x, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jpp = zeros(m,o,o,n);
    dfdpdphan = data.dfdpdphan;
    for i=1:n
      Jpp(:,:,:,i) = dfdpdphan(x(:,i), p(:,i));
    end
  end
end
end

function Jtx = aut_DFDTDX(~, ~, x, ~)
%DFDTDX   Vectorized evaluation of autonomous d2F/dtdx.
Jtx = zeros([size(x,1), size(x)]);
end

function Jtp = aut_DFDTDP(~, ~, x, p)
%DFDTDX   Vectorized evaluation of autonomous d2F/dtdp.
Jtp = zeros([size(x,1), size(p)]);
end

function Jtt = aut_DFDTDT(~, ~, x, ~)
%DFDTDT   Vectorized evaluation of autonomous d2F/dtdt.
Jtt = zeros(size(x));
end

function f = het_F(data, t, x, p)
%F   Vectorized evaluation of non-autonomous F.

if data.ode.vectorized
  f = data.fhan(t, x, p);
else
  [m, n] = size(x);
  f = zeros(m,n);
  fhan = data.fhan;
  for i=1:n
    f(:,i) = fhan(t(i), x(:,i), p(:,i));
  end
end
end

function Jx = het_DFDX(data, t, x, p)
%DFDX   Vectorized evaluation of non-autonomous dF/dx.

if data.ode.vectorized
  if isempty(data.dfdxhan)
    f  = @(x,p) data.fhan(p(1,:), x, p(2:end,:));
    tp = [ t ; p ];
    Jx = coco_ezDFDX('f(x,p)v', f, x, tp);
  else
    Jx = data.dfdxhan(t, x, p);
  end
else
  if isempty(data.dfdxhan)
    f  = @(x,p) data.fhan(p(1,:), x, p(2:end,:));
    tp = [ t ; p ];
    Jx = coco_ezDFDX('f(x,p)', f, x, tp);
  else
    [m, n] = size(x);
    Jx = zeros(m,m,n);
    dfdxhan = data.dfdxhan;
    for i=1:n
      Jx(:,:,i) = dfdxhan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jp = het_DFDP(data, t, x, p)
%DFDP   Vectorized evaluation of non-autonomous dF/dp.

if data.ode.vectorized
  if isempty(data.dfdphan)
    f  = @(x,p) data.fhan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jp = coco_ezDFDP('f(x,p)v', f, tx, p);
  else
    Jp = data.dfdphan(t, x, p);
  end
else
  if isempty(data.dfdphan)
    f  = @(x,p) data.fhan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jp = coco_ezDFDP('f(x,p)', f, tx, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jp = zeros(m,o,n);
    dfdphan = data.dfdphan;
    for i=1:n
      Jp(:,:,i) = dfdphan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jt = het_DFDT(data, t, x, p)
%DFDT   Vectorized evaluation of non-autonomous dF/dt.

if data.ode.vectorized
  if isempty(data.dfdthan)
    [m, n] = size(x);
    f  = @(x,p) data.fhan(x(1,:), p(1:m,:), p(m+1:end,:));
    xp = [ x ; p ];
    Jt = reshape(coco_ezDFDX('f(x,p)v', f, t, xp), [m n]);
  else
    Jt = data.dfdthan(t, x, p);
  end
else
  if isempty(data.dfdthan)
    [m, n] = size(x);
    f  = @(x,p) data.fhan(x(1,:), p(1:m,:), p(m+1:end,:));
    xp = [ x ; p ];
    Jt = reshape(coco_ezDFDX('f(x,p)', f, t, xp), [m n]);
  else
    [m, n] = size(x);
    Jt = zeros(m,n);
    dfdthan = data.dfdthan;
    for i=1:n
      Jt(:,i) = dfdthan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jxx = het_DFDXDX(data, t, x, p)
%DFDXDX   Vectorized evaluation of non-autonomous d2F/dx2.

if data.ode.vectorized
  if isempty(data.dfdxdxhan)
    f  = @(x,p) data.dfdxhan(p(1,:), x, p(2:end,:));
    tp = [ t ; p ];
    Jxx = coco_ezDFDX('f(x,p)v', f, x, tp);
  else
    Jxx = data.dfdxdxhan(t, x, p);
  end
else
  if isempty(data.dfdxdxhan)
    f  = @(x,p) data.dfdxhan(p(1,:), x, p(2:end,:));
    tp = [ t ; p ];
    Jxx = coco_ezDFDX('f(x,p)', f, x, tp);
  else
    [m, n] = size(x);
    Jxx = zeros(m,m,m,n);
    dfdxdxhan = data.dfdxdxhan;
    for i=1:n
      Jxx(:,:,:,i) = dfdxdxhan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jxp = het_DFDXDP(data, t, x, p)
%DFDXDP   Vectorized evaluation of non-autonomous d2F/dxdp.

if data.ode.vectorized
  if isempty(data.dfdxdphan)
    f  = @(x,p) data.dfdxhan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jxp = coco_ezDFDP('f(x,p)v', f, tx, p);
  else
    Jxp = data.dfdxdphan(t, x, p);
  end
else
  if isempty(data.dfdxdphan)
    f  = @(x,p) data.dfdxhan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jxp = coco_ezDFDP('f(x,p)', f, tx, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jxp = zeros(m,m,o,n);
    dfdxdphan = data.dfdxdphan;
    for i=1:n
      Jxp(:,:,:,i) = dfdxdphan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jpp = het_DFDPDP(data, t, x, p)
%DFDXDP   Vectorized evaluation of non-autonomous d2F/dxdp.

if data.ode.vectorized
  if isempty(data.dfdpdphan)
    f  = @(x,p) data.dfdphan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jpp = coco_ezDFDP('f(x,p)v', f, tx, p);
  else
    Jpp = data.dfdpdphan(t, x, p);
  end
else
  if isempty(data.dfdpdphan)
    f  = @(x,p) data.dfdphan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jpp = coco_ezDFDP('f(x,p)', f, tx, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jpp = zeros(m,o,o,n);
    dfdpdphan = data.dfdpdphan;
    for i=1:n
      Jpp(:,:,:,i) = dfdpdphan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jtx = het_DFDTDX(data, t, x, p)
%DFDTDX   Vectorized evaluation of non-autonomous d2F/dtdx.

if data.ode.vectorized
  if isempty(data.dfdtdxhan)
    f  = @(x,p) data.dfdthan(p(1,:), x, p(2:end,:));
    tp = [ t ; p ];
    Jtx = coco_ezDFDX('f(x,p)v', f, x, tp);
  else
    Jtx = data.dfdtdxhan(t, x, p);
  end
else
  if isempty(data.dfdtdxhan)
    f  = @(x,p) data.dfdthan(p(1,:), x, p(2:end,:));
    tp = [ t ; p ];
    Jtx = coco_ezDFDX('f(x,p)', f, x, tp);
  else
    [m, n] = size(x);
    Jtx = zeros(m,m,n);
    dfdtdxhan = data.dfdtdxhan;
    for i=1:n
      Jtx(:,:,i) = dfdtdxhan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jtp = het_DFDTDP(data, t, x, p)
%DFDTDP   Vectorized evaluation of non-autonomous d2F/dtdp.

if data.ode.vectorized
  if isempty(data.dfdtdphan)
    f  = @(x,p) data.dfdthan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jtp = coco_ezDFDP('f(x,p)v', f, tx, p);
  else
    Jtp = data.dfdtdphan(t, x, p);
  end
else
  if isempty(data.dfdtdphan)
    f  = @(x,p) data.dfdthan(x(1,:), x(2:end,:), p);
    tx = [ t ; x ];
    Jtp = coco_ezDFDP('f(x,p)', f, tx, p);
  else
    [m, n] = size(x);
    o  = size(p,1);
    Jtp = zeros(m,o,n);
    dfdtdphan = data.dfdtdphan;
    for i=1:n
      Jtp(:,:,i) = dfdtdphan(t(i), x(:,i), p(:,i));
    end
  end
end
end

function Jtt = het_DFDTDT(data, t, x, p)
%DFDTDT   Vectorized evaluation of non-autonomous d2F/dt2.

if data.ode.vectorized
  if isempty(data.dfdtdthan)
    [m, n] = size(x);
    f  = @(x,p) data.dfdthan(x(1,:), p(1:m,:), p(m+1:end,:));
    xp = [ x ; p ];
    Jtt = reshape(coco_ezDFDX('f(x,p)v', f, t, xp), [m n]);
  else
    Jtt = data.dfdtdthan(t, x, p);
  end
else
  if isempty(data.dfdtdthan)
    [m, n] = size(x);
    f  = @(x,p) data.dfdthan(x(1,:), p(1:m,:), p(m+1:end,:));
    xp = [ x ; p ];
    Jtt = reshape(coco_ezDFDX('f(x,p)', f, t, xp), [m n]);
  else
    [m, n] = size(x);
    Jtt = zeros(m,m,n);
    dfdtdthan = data.dfdtdthan;
    for i=1:n
      Jtt(:,:,i) = dfdtdthan(t(i), x(:,i), p(:,i));
    end
  end
end
end
