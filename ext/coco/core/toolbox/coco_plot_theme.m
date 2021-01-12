function theme = coco_plot_theme(tb_info)
%COCO_PLOT_THEME    Default plot theme
%
%   COCO_PLOT_THEME with no arguments returns a struct with global default
%   settings for all toolboxes.
%
% See also: COCO_PLOT_BD, COCO_PLOT_SOL

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_plot_theme.m 3015 2017-08-19 03:25:33Z hdankowicz $

theme = struct();
theme.bd.col1  = '';
theme.bd.col2  = '';
theme.sol.col1 = '';
theme.sol.col2 = '';
theme.lspec    = {'b-', 'LineWidth', 1};
theme.ustab    = '';
theme.usept    = {};
theme.special  = {};
theme.FP       = {'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 8};
theme.BP       = {'kx', 'LineWidth', 2, 'Color', 'b', 'MarkerSize', 8};
theme.EP       = {'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5};
theme.axlabs   = true;
theme.xlab     = '';
theme.ylab     = '';
theme.zlab     = '';
theme.xlstyle  = {'FontSize', 12};
theme.ylstyle  = {'FontSize', 12};
theme.zlstyle  = {'FontSize', 12};

theme.plot_sol = @default_plot_sol;

if nargin>=1
  TB = tb_info.tb;
  thmfunc = str2func(sprintf('%s_plot_theme', TB));
  thm = thmfunc(tb_info.(TB).branch_type);
  theme = coco_merge(theme, thm);
end
end

function thm = default_plot_sol(varargin)
error('%s: no plot function defined for plotting solutions\n');
end
