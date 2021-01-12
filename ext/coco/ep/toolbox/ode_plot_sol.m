function ode_plot_sol(varargin)
%ODE_PLOT_SOL    Plot solutions.
%
%   ODE_PLOT_SOL([THEME], RUN, [LABS], [COL1, [IDX1], [COL2, [IDX2], [COL3, [IDX3]]]])
%
%   ODE_PLOT_SOL(RUN) plots the first against the second dimension of all
%   labelled solutions computed along branch with name RUN.
%
%   ODE_PLOT_SOL(RUN, LABS) plots the solutions with labels LABS.
%
%   ODE_PLOT_SOL(..., COL1, [IDX1], ...) plots the specified dimension of
%   bifurcation data or solution components. Columns can be:
%     'x'  - solution vector
%     't'  - time (PO only)
%     'tn' - normalised time (PO only)
%   or any column from the bifurcation data. For columns with
%   multi-dimensional data IDX1 ... IDX3 select the dimension to plot. A
%   column might be specified several times. ODE_PLOT_SOL tries to be smart
%   when selecting default dimensions to plot.
%
%   ODE_PLOT_SOL(THEME, ...) uses a user-specified theme for plotting. THEME
%   is a struct and the fields given in THEME are merged into the default
%   theme recursively. See ode_plot_theme, ep_plot_theme and po_plot_theme
%   for help on default themes and what fields are available.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_plot_sol.m 2946 2016-04-07 12:54:52Z fschild $

thm = struct();

grammar   = '[THM] RUN [LABS] [C1 [I1] [C2 [I2] [C3 [I3]]]]';
args_spec = {
  'THM',     '', 'stct',  'theme',    thm, 'read', {}
  'RUN',     '',  'str',    'run',     '', 'read', {}
  'LABS',    '',  'num',   'labs',  'ALL', 'read', {}
   'C1',     '',  'str',   'col1',    'x', 'read', {}
   'I1',     '',  'num',   'idx1',     [], 'read', {}
   'C2',     '',  'str',   'col2',    'x', 'read', {}
   'I2',     '',  'num',   'idx2',     [], 'read', {}
   'C3',     '',  'str',   'col3',     '', 'read', {}
   'I3',     '',  'num',   'idx3',     [], 'read', {}
  };
opts_spec = {};
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:}); %#ok<NASGU>

[bd tb_info] = coco_bd_read(args.run, 'bd', 'tb_info', 'bddat');
thm = coco_merge(ode_plot_theme(tb_info), args.theme);

if isempty(args.col1)
  args.col1 = 'x';
end

if isempty(args.idx1)
  args.idx1 = 1;
end

if isempty(args.col2)
  args.col2 = 'x';
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
else
  dim = 3;
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

if isempty(thm.xlab)
  thm.xlab = args.col1;
end
if isempty(thm.ylab)
  thm.ylab = args.col2;
end
if isempty(thm.zlab)
  thm.zlab = args.col3;
end

if ischar(args.labs) || isempty(args.labs)
  args.labs = coco_bd_labs(bd, 'all');
end

HH = ~ishold;

thm.plot_sol(thm, args, tb_info, dim, bd, HH);

if thm.axlabs
  xlabel(thm.xlab, thm.xlstyle{:});
  ylabel(thm.ylab, thm.ylstyle{:});
  if dim==3
    zlabel(thm.zlab, thm.zlstyle{:});
  end
end

if HH; hold off; end

end
