function coco_plot_sol(varargin)
%COCO_PLOT_SOL    Plot solutions.
%
%   COCO_PLOT_SOL([THEME], RUN, [LABS], OID, [OIDX], [COL1, [IDX1], [COL2, [IDX2], [COL3, [IDX3]]]])
%
%   COCO_PLOT_SOL(RUN) plots the first against the second dimension of all
%   labelled solutions computed along branch with name RUN.
%
%   COCO_PLOT_SOL(RUN, LABS) plots the solutions with labels LABS.
%
%   COCO_PLOT_SOL(..., COL1, [IDX1], ...) plots the specified dimension of
%   bifurcation data or solution components. Columns can be toolbox
%   specific or any column from the bifurcation data. For columns with
%   multi-dimensional data IDX1 ... IDX3 select the dimension to plot. A
%   column might be specified several times. COCO_PLOT_SOL tries to be
%   smart when selecting default dimensions to plot.
%
%   COCO_PLOT_SOL(THEME, ...) uses a user-specified theme for plotting.
%   THEME is a struct and the fields given in THEME are merged into the
%   default theme recursively. See coco_plot_theme for help on default
%   themes and what fields are available.
%
% See also: COCO_PLOT_THEME

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_plot_sol.m 2946 2016-04-07 12:54:52Z fschild $

thm = struct();

grammar   = '[THM] RUN [LABS] OID [OIDX] [COL1 [IDX1] [COL2 [IDX2] [COL3 [IDX3]]]]';
args_spec = {
   'THM',     '',   'stct', 'theme',   thm, 'read', {}
   'RUN',     '',    'str',   'run',    '', 'read', {}
  'LABS',     '',  '[num]',  'labs', 'ALL', 'read', {}
   'OID',     '',    'str',   'oid',    '', 'read', {}
  'OIDX',     '',  '[num]',  'oidx',    [], 'read', {}
  'COL1', 'cell',  '{str}',  'col1',    '', 'read', {}
  'IDX1',     '', 'num|@f',  'idx1',    [], 'read', {}
  'COL2', 'cell',  '{str}',  'col2',    '', 'read', {}
  'IDX2',     '', 'num|@f',  'idx2',    [], 'read', {}
  'COL3', 'cell',  '{str}',  'col3',    '', 'read', {}
  'IDX3',     '', 'num|@f',  'idx3',    [], 'read', {}
  };
opts_spec = {};
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});  %#ok<ASGLU>

if ~isempty(args.oidx)
  oid = cell(1,numel(args.oidx));
  for i=1:numel(args.oidx)
    oid{i} = sprintf('%s%d', args.oid, args.oidx(i));
  end
  args.oid = oid;
else
  args.oid = {args.oid};
end
tb_info = coco_bd_read(args.run, coco_get_id(args.oid{1}, 'tb_info'));
if ~isempty(tb_info)
  thm = coco_plot_theme(tb_info);
else
  thm = coco_plot_theme();
end
thm = coco_merge(thm, args.theme);

if isempty(args.col1)
  args.col1 = {thm.sol.col1};
end

if isempty(args.idx1)
  args.idx1 = 1;
end

if isempty(args.col2)
  args.col2 = {thm.sol.col2};
end

if isempty(args.idx2)
  args.idx2 = 1;
  if numel(args.col1)==1 && isnumeric(args.idx1) ...
      && strcmp(args.col2, args.col1)
    args.idx2 = args.idx1 + 1;  
  end
end

if isempty(args.col3)
  dim = 2;
else
  dim = 3;
  if isempty(args.idx3)
    args.idx3 = 1;
    if numel(args.col1)==1 && numel(args.col2)==1 ...
        && isnumeric(args.idx1) && isnumeric(args.idx2)
      flags = strcmp(args.col3, {args.col1, args.col2});
      if all(flags)
        args.idx3 = max(args.idx1, args.idx2)+1;
      elseif flags(1)
        args.idx3 = args.idx1+1;
      elseif flags(2)
        args.idx3 = args.idx2+1;
      end
    elseif numel(args.col1)==1 && isnumeric(args.idx1) && ...
        strcmp(args.col3, args.col1)
      args.idx3 = args.idx1 + 1;
    elseif numel(args.col2)==1 && isnumeric(args.idx2) && ...
        strcmp(args.col3, args.col2)
      args.idx3 = args.idx2 + 1;
    end
  end
end

if ~isa(args.idx1, 'function_handle')
  assert(numel(args.col1)==1, '%s: %s', ...
    'incorrect indexing along 1st dimension', mfilename);
end

if ~isa(args.idx2, 'function_handle')
  assert(numel(args.col2)==1, '%s: %s', ...
    'incorrect indexing along 2nd dimension', mfilename);
end

if dim==3 && ~isa(args.idx3, 'function_handle')
  assert(numel(args.col3)==1, '%s: %s', ...
    'incorrect indexing along 3rd dimension', mfilename);
end

bd = coco_bd_read(args.run, 'bd');
if ischar(args.labs) || isempty(args.labs)
  args.labs = coco_bd_labs(bd, 'all');
end

HH = ~ishold;

thm = thm.plot_sol(thm, args, dim, bd, HH);

if thm.axlabs
  xlabel(thm.xlab, thm.xlstyle{:});
  ylabel(thm.ylab, thm.ylstyle{:});
  zlabel(thm.zlab, thm.zlstyle{:});
end

if HH; hold off; end

end
