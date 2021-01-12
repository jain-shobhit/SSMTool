function ode_plot_bd(varargin)
%ODE_PLOT_BD    Plot bifurcation diagram.
%
%   ODE_PLOT_BD([THEME], RUN, [COL1, [IDX1], [COL2, [IDX2], [COL3, [IDX3]]]])
%  
%   ODE_PLOT_BD(RUN) plots the default solution measure ||x|| as a function
%   of the primary continuation parameter. RUN is the run name of the
%   branch.
%
%   ODE_PLOT_BD(RUN, COL1, [IDX1], ...) plots the user specified columns of
%   the bifurcation data. For columns with multi-dimensional data IDX1 ...
%   IDX3 select the dimension to plot. A column might be specified several
%   times. ODE_PLOT_BD tries to be smart when selecting default dimensions
%   to plot.
%
%   ODE_PLOT_BD(THEME, ...) uses a user-specified theme for plotting. THEME
%   is a struct and the fields given in THEME are merged into the default
%   theme recursively. See ode_plot_theme, ep_plot_theme and po_plot_theme
%   for help on default themes and what fields are available.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_plot_bd.m 3015 2017-08-19 03:25:33Z hdankowicz $

thm = struct();

grammar   = '[THM] RUN [C1 [I1] [C2 [I2] [C3 [I3]]]]';
args_spec = {
  'THM',     '', 'stct',  'theme', thm, 'read', {}
  'RUN',     '',  'str',    'run',  '', 'read', {}
   'C1',     '',  'str',   'col1',  '', 'read', {}
   'I1',     '',  'num',   'idx1',  [], 'read', {}
   'C2',     '',  'str',   'col2',  '', 'read', {}
   'I2',     '',  'num',   'idx2',  [], 'read', {}
   'C3',     '',  'str',   'col3',  '', 'read', {}
   'I3',     '',  'num',   'idx3',  [], 'read', {}
  };
opts_spec = {};
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:}); %#ok<NASGU>

[bd, tb_info, bddat] = coco_bd_read(args.run, 'bd', 'tb_info', 'bddat');
thm = ode_plot_theme(tb_info);
if isfield(args.theme, 'special')
  args.theme.special = union(thm.special, args.theme.special);
else
  args.theme.special = thm.special;
end
thm = coco_merge(thm, args.theme);

if isempty(args.col1)
  if isempty(thm.bd.col1)
    args.col1 = bddat.op_names{1};
  else
    args.col1 = thm.bd.col1;
  end
end

if isempty(args.idx1)
  args.idx1 = 1;
end

if isempty(args.col2)
  if isempty(thm.bd.col2)
    error('%s: theme error: ''theme.bd.col2'' empty');
  else
    args.col2 = thm.bd.col2;
  end
end

if isempty(args.idx2)
  if strcmp(args.col1, args.col2)
    args.idx2 = args.idx1+1;
  else
    args.idx2 = 1;
  end
end

if isempty(args.col3)
  dim = 2;
  spfunc = @(HH,S,PT,SP, x,y,z, ST1,ST2) stab_plot(HH,S,PT,SP, x,y, ST1,ST2);
  pfunc  = @(x,y,z,varargin) plot(x,y, varargin{:});
else
  dim = 3;
  spfunc = @(HH,S,PT,SP, x,y,z, ST1,ST2) stab_plot3(HH,S,PT,SP, x,y,z, ST1,ST2);
  pfunc  = @(x,y,z,varargin) plot3(x,y,z, varargin{:});
  if isempty(args.idx3)
    flags = strcmp(args.col3, {args.col1, args.col2});
    if all(flags)
      args.idx3 = max(args.idx1, args.idx2)+1;
    elseif flags(1)
      args.idx3 = args.idx1+1;
    elseif flags(2)
      args.idx3 = args.idx2+1;
    else
      args.idx3 = 1;
    end
  end
end

x = coco_bd_col(bd, args.col1);
x = x(args.idx1,:);
if isempty(thm.xlab)
  thm.xlab = args.col1;
end
y = coco_bd_col(bd, args.col2);
y = y(args.idx2,:);
if isempty(thm.ylab)
  thm.ylab = args.col2;
end
if dim == 2
  z = nan(size(x));
else
  z = coco_bd_col(bd, args.col3);
  z = z(args.idx3,:);
  if isempty(thm.zlab)
    thm.zlab = args.col3;
  end
end

if isempty(thm.ustab)
  S   = zeros(size(x));
  PT  = cell(size(x));
  SP  = {};
  ST1 = thm.lspec;
  ST2 = thm.lspec;
else
  S   = coco_bd_col(bd, thm.ustab);
  PT  = coco_bd_col(bd, 'TYPE', false);
  SP  = thm.usept;
  ST1 = thm.lspec_s;
  ST2 = thm.lspec_u;
end

HH = ~ishold;

spfunc(HH,S,PT,SP, x,y,z, ST1,ST2);

N = numel(thm.special);
for i=1:N
  SP = thm.special{i};
  LS = thm.(SP);
  idx = coco_bd_idxs(bd, SP);
  pfunc(x(idx),y(idx),z(idx), LS{:});
end

if thm.axlabs
  xlabel(thm.xlab, thm.xlstyle{:});
  ylabel(thm.ylab, thm.ylstyle{:});
  if dim==3
    zlabel(thm.zlab, thm.zlstyle{:});
  end
end

if HH; hold off; end

end

function stab_plot(HH,S,PT,SP, x,y, ST1,ST2)
CF  = @(x) (isempty(x) || any(strcmp(x, {'RO', 'EP'})));
PTF = ~cellfun(CF, PT);
NSP = numel(SP);
II  = 1;
EI  = numel(x);
I   = II;
SS  = (S(II)==0);
while true
  if II>=EI; break; end
  [I I0] = next_index(EI, SS, S, PT, PTF, NSP, SP, I, II);
  if SS
    plot(x(II:I0), y(II:I0), ST1{:});
  else
    plot(x(II:I0), y(II:I0), ST2{:});
  end
  II=I0;
  SS=(S(I)==0);
  if HH; hold on; end
end

end

function stab_plot3(HH,S,PT,SP, x,y,z, ST1,ST2)
CF  = @(x) (isempty(x) || any(strcmp(x, {'RO', 'EP'})));
PTF = ~cellfun(CF, PT);
NSP = numel(SP);
II = 1;
EI = numel(x);
I  = II;
SS = (S(II)==0);
while true
  if II>=EI; break; end
  [I I0] = next_index(EI, SS, S, PT, PTF, NSP, SP, I, II);
  if SS
    plot3(x(II:I0), y(II:I0), z(II:I0), ST1{:});
  else
    plot3(x(II:I0), y(II:I0), z(II:I0), ST2{:});
  end
  II=I0;
  SS=(S(I)==0);
  if HH; hold on; end
end

end

function [I I0] = next_index(EI, SS, S, PT, PTF, NSP, SP, I, II)
while (I<EI && SS==(S(I)==0)); I=I+1; end
I0 = I-1;
while (I0>II && PTF(I0)); I0 = I0-1; end
eflag = false;
for J=1:NSP
  CSP = SP{J};
  for K=(I0+1):(I-1)
    if strcmp(PT{K},CSP)
      I0 = K;
      eflag = true;
      break;
    end
  end
  if eflag; break; end
end
if ~eflag
  I0 = I;
end
end
