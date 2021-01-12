function prob = coco_change_func(prob, data, varargin)
% COCO_CHANGE_FUNC   COCO function reconstructor following remeshing
%
% PROB = COCO_CHANGE_FUNC(PROB, DATA, VARARGIN)
% VARARGIN = OPTS...
%
% OPTS = { ('uidx'|'xidx') ('all'|I) }
% OPTS = { ('u0'|'x0') U0 }
% OPTS = { ('fdim'|'ydim') N }
% OPTS = { 'vecs' V0 }

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coco_change_func.m 3100 2019-06-07 17:31:58Z hdankowicz $

efunc = prob.efunc;
func  = efunc.funcs(prob.efunc.cfidx);
type  = func.type;

switch type
  case 'zero'
    pnum = 0;
  otherwise
    pnum   = func.pnum;
    pnames = func.pnames;
end

efopts.x0    = [];
efopts.xidx  = [];
efopts.V0    = [];
efopts.props = func.props;
xname        = 'u0';

argidx = 1;
while argidx+2<=nargin
  oarg   = varargin{argidx};
  argidx = argidx + 1;
  oname  = lower(oarg);
  switch oname
      
    case { 'uidx' 'xidx' }
      efopts.xidx = varargin{argidx};
      efopts.xidx = efopts.xidx(:)';
      argidx      = argidx + 1;
    
    case { 'u0' 'x0' }
      xname       = oarg;
      efopts.x0   = varargin{argidx};
      argidx      = argidx + 1;
      
    case { 'fdim' 'ydim' }
      efopts.fdim = varargin{argidx};
      argidx      = argidx + 1;
      
    case 'vecs'
      efopts.V0   = varargin{argidx};
      argidx      = argidx + 1;
      
    case { 'prop' 'property' }
      name        = varargin{argidx};
      val         = varargin{argidx+1};
      argidx      = argidx + 2;
      efopts.props.(name) = val;
      
    otherwise
      if ischar(oarg)
        error('%s: option ''%s'' not recognised', mfilename, oarg);
      else
        error('%s: in argument %d: expected string, got a ''%s''', ...
          mfilename, argidx-1, class(oarg));
      end
  end
end

xnum = numel(efopts.x0);
if numel(efopts.x0)~=size(efopts.V0,1)
  emsg = sprintf('%s:', mfilename);
  emsg = sprintf('%s dimension of vectors in vecs [%d] does not', ...
    emsg, size(efopts.V0,1));
  emsg = sprintf('%s\nagree with dimension of %s [%d]', ...
    emsg, xname, xnum);
  error(emsg); %#ok<SPERR>
end

if ~strcmpi(efopts.xidx, 'all')
  x0          = efopts.x0(:);
  xidx        = efunc.x_dim + (1:xnum);
  efunc.x_dim = efunc.x_dim + xnum;
  efunc.x_idx = [ efunc.x_idx             xidx ];
  efopts.xidx = [ efopts.xidx             xidx ];
  efunc.x0    = [ efunc.x0 ;              x0   ];
  x0          = [ efunc.x0(efopts.xidx) ; x0   ];
  efunc.V0    = [ efunc.V0 ;         efopts.V0 ];
end

%% create function object and update efunc.f_dim

func.data  = data;
func.x_idx = efopts.xidx;
func.props = efopts.props;

if ~isfield(efopts, 'fdim')
  switch type
    case 'zero'
      if ~strcmpi(efopts.xidx, 'all')
        [prob, func.data, efunc.chart, f] = ...
          efunc_call_F(prob, func.data, efunc.chart, func, x0);
        efopts.fdim = numel(f);
      else
        efopts.fdim = numel(func.f_idx);
      end
    otherwise
      efopts.fdim = pnum;
  end
end

fdim = efopts.fdim;
switch type
  case { 'zero' 'inactive' 'active' 'internal' 'inequality' }
    fidx             = efunc.f_dim + (1:fdim);
    efunc.f_dim      = efunc.f_dim + fdim;
    efunc.pidx2fidx  = [ efunc.pidx2fidx ; fidx' ];
  otherwise
    fidx            = [];
    efunc.pidx2fidx = [ efunc.pidx2fidx ; zeros(fdim,1) ];
end
func.f_idx  = fidx;
efunc.f_idx = [efunc.f_idx fidx];

switch type
  case 'zero'
    midx            = [];
    pnames          = cell(1,fdim);
    efunc.pidx2midx = [ efunc.pidx2midx ; zeros(fdim,1) ];
  otherwise
    midx            = efunc.m_dim + (1:fdim);
    efunc.m_dim     = efunc.m_dim + fdim;
    efunc.pidx2midx = [ efunc.pidx2midx ; midx' ];
end
func.m_idx          = midx;

func.req_idx = midx;
for i=1:numel(func.requires)
  idx = strcmpi(func.requires{i}, efunc.identifyers);
  efunc.funcs(idx).req_idx = [efunc.funcs(idx).req_idx midx];
end

efunc.funcs(prob.efunc.cfidx) = func;

%% add pnames to list of parameters and update efunc.p_dim

idx              = numel(efunc.idx2par) + (1:pnum);
efunc.idx2par    = [ efunc.idx2par pnames ];
pararray         = sprintf('%s_pars', type);
efunc.(pararray) = [ efunc.(pararray) idx ];

switch type
	case { 'active' 'internal' 'inequality' }
    efunc.p_dim = efunc.p_dim + pnum;
end

%% return efunc in prob

prob.efunc = efunc;

end
