function prob = coco_change_comp(prob, data, varargin)
%COCO_CHANGE_COMP   Complementary function reconstructor following remeshing
%
% PROB = COCO_CHANGE_COMP(PROB, DATA, VARARGIN)
% VARARGIN = OPTS...
%
% OPTS = { 'uidx' I }
% OPTS = { 'lidx' I }
% OPTS = { 'vidx' I }
% OPTS = { 'v0' v0 }
% OPTS = { 'vdim' N }
% OPTS = { 'vecs' V0 }

% Copyright (C) Harry Dankowicz, Mingwu Li

efunc         = prob.efunc;
adjoint       = prob.adjoint;
complementary = prob.complementary;

cfopts.uidx = [];
cfopts.lidx = [];
cfopts.vidx = [];
cfopts.v0   = [];
cfopts.V0   = [];

argidx = 1;
while(argidx<=nargin-2)
  oarg   = varargin{argidx};
  argidx = argidx + 1;
  oname  = lower(oarg);
  switch oname
    
    case 'v0'
      cfopts.v0 = varargin{argidx};
      argidx = argidx + 1;

    case 'uidx'
      cfopts.uidx = varargin{argidx};
      cfopts.uidx = cfopts.uidx(:)';
      argidx = argidx + 1;
      
    case 'lidx'
      cfopts.lidx = varargin{argidx};
      cfopts.lidx = cfopts.lidx(:)';
      argidx = argidx + 1;
      
    case 'vidx'
      cfopts.vidx = varargin{argidx};
      cfopts.vidx = cfopts.vidx(:)';
      argidx = argidx + 1;

    case 'vdim'
      cfopts.vdim = varargin{argidx};
      argidx = argidx + 1;

    case 'vecs'
      cfopts.V0 = varargin{argidx};
      argidx = argidx + 1;

    otherwise
      if ischar(oarg)
        error('%s: option ''%s'' not recognised', mfilename, oarg);
      else
        error('%s: in argument %d: expected string, got a ''%s''', ...
          mfilename, 1+argidx, class(oarg));
      end
  end
  
end

%% create function object and update complementary function

u0 = efunc.x0(cfopts.uidx);
l0 = adjoint.l0(cfopts.lidx);
v0 = cfopts.v0(:);
vnum = numel(cfopts.v0);

if vnum>0
  vidx = complementary.v_dim + (1:vnum);
else
  vidx = [];
end
complementary.v_dim = complementary.v_dim + vnum;
complementary.v0    = [ complementary.v0              ;   v0 ];
v0                  = [ complementary.v0(cfopts.vidx) ;   v0 ];
complementary.v_idx = [ complementary.v_idx             vidx ];
cfopts.vidx         = [ cfopts.vidx                     vidx ];

%% create function object and update adjoint

comp = complementary.funcs(prob.complementary.cfidx);

comp.u_idx = cfopts.uidx;
comp.l_idx = cfopts.lidx;
comp.v_idx = cfopts.vidx;
comp.v_num = vnum;

if ~isfield(cfopts, 'fdim')
  [data, f] = complementary_call_F(prob, data, comp, u0, l0, v0);
  cfopts.fdim = numel(f);
end
comp.data = data;

fdim = cfopts.fdim;
switch comp.type
  case { 'zero' }
    fidx                = complementary.f_dim + (1:fdim);
    complementary.f_dim = complementary.f_dim + fdim;
  case { 'inactive' 'active' 'internal' }
    fidx                = complementary.f_dim + (1:fdim);
    complementary.f_dim = complementary.f_dim + fdim;
  otherwise
    fidx                = [];
end
comp.f_idx = fidx;
complementary.f_idx = [complementary.f_idx fidx];

complementary.funcs(prob.complementary.cfidx) = comp;

complementary.V0 = [ complementary.V0 ; cfopts.V0 ];

prob.complementary = complementary;

end
