function prob = coco_change_adjt(prob, data, varargin)
%COCO_CHANGE_ADJT   Adjoint reconstructor following remeshing
%
% PROB = COCO_CHANGE_ADJT(PROB, DATA, VARARGIN)
% VARARGIN = OPTS...
%
% OPTS = { 'l0' L0 }
% OPTS = { 'aidx' I }
% OPTS = { 'adim' [N,N] }
% OPTS = { 'vecs' VL0 }

% Copyright (C) Harry Dankowicz, Mingwu Li

adjoint = prob.adjoint;
adjt    = adjoint.funcs(prob.adjoint.cfidx);
fid     = adjt.identifyer;

efunc = prob.efunc;
idx   = find(strcmpi(fid, efunc.identifyers), 1);
x0    = coco_get_func_data(prob, fid, 'x0');
func  = efunc.funcs(idx);
x_idx = func.x_idx;

afopts.aidx = [];
afopts.l0   = [];
afopts.Vl0  = [];

argidx = 1;
while(argidx<=nargin-2)
  oarg   = varargin{argidx};
  argidx = argidx + 1;
  oname  = lower(oarg);
  switch oname
    
    case 'l0'
      afopts.l0 = varargin{argidx};
      argidx = argidx + 1;

    case 'aidx'
      afopts.aidx = varargin{argidx};
      afopts.aidx = afopts.aidx(:)';
      argidx = argidx + 1;

    case 'adim'
      afopts.adim = varargin{argidx};
      argidx = argidx + 1;

    case 'vecs'
      afopts.Vl0 = varargin{argidx};
      argidx     = argidx + 1;

    otherwise
      if ischar(oarg)
        error('%s: option ''%s'' not recognised', mfilename, oarg);
      else
        error('%s: in argument %d: expected string, got a ''%s''', ...
          mfilename, 1+argidx, class(oarg));
      end
  end
  
end

%% create function object and update adjoint

if ~isfield(afopts, 'adim')
  [data, J] = adjt.F(prob, data, x0);
  afopts.adim = size(J);
end
afdim = afopts.adim(1);
axdim = afopts.adim(2);
axnum = numel(afopts.aidx);

afidx          = adjoint.a_dim(1) + (1:afdim);
axidx          = adjoint.a_dim(2) + (1:axdim-axnum);
adjoint.a_dim  = adjoint.a_dim + [afdim, axdim-axnum];
adjoint.af_idx = [ adjoint.af_idx afidx ];
adjoint.ax_idx = [ adjoint.ax_idx axidx ];
axidx          = [ afopts.aidx    axidx ];

adjt.data       = data;
adjt.af_idx     = afidx;
adjt.ax_idx     = axidx;
adjt.x_idx      = x_idx;

adjoint.funcs(prob.adjoint.cfidx) = adjt;

adjoint.l0  = [ adjoint.l0  ;  afopts.l0 ];
adjoint.Vl0 = [ adjoint.Vl0 ; afopts.Vl0 ];
switch func.type
  case 'zero'
  case 'inequality'
    adjoint.s_idx  = [ adjoint.s_idx ; adjt.af_idx' ];
  otherwise
    adjoint.l_idx  = [ adjoint.l_idx ; adjt.af_idx' ];
end

prob.adjoint = adjoint;

end
