function theme = ode_plot_theme(tb_info)
%ODE_PLOT_THEME    Default plot theme for ODE_PLOT functions.
%
%   ODE_PLOT_THEME with no arguments returns a struct with global default
%   settings for all ODE toolboxes.
%
%   ODE_PLOT_THEME(TB_INFO) computes a default theme struct by merging the
%   global ODE defaults with toolbox specific defaults. The toolbox is
%   identified in the toolbox info struct TB_INFO that every toolbox in the
%   ode toolbox family attaches to the bifurcation data.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_plot_theme.m 3015 2017-08-19 03:25:33Z hdankowicz $

theme = struct();
theme.bd.col1  = '';
theme.bd.col2  = '';
theme.sol.col1 = '';
theme.sol.col2 = '';
theme.lspec    = {'b-', 'LineWidth', 1};
theme.lspec_s  = theme.lspec;
theme.lspec_u  = {'b--', 'LineWidth', 1};
theme.ustab    = '';
theme.usept    = {};
theme.special  = {'FP', 'BP'};
theme.FP       = {'kd', 'MarkerFaceColor', 'w', 'MarkerSize', 8};
theme.BP       = {'kx', 'LineWidth', 2, 'Color', 'b', 'MarkerSize', 8};
theme.axlabs   = true;
theme.xlab     = '';
theme.ylab     = '';
theme.zlab     = '';
theme.xlstyle  = {'FontSize', 12};
theme.ylstyle  = {'FontSize', 12};
theme.zlstyle  = {'FontSize', 12};

c = [1 1 0.8];
theme.labs   = false;
theme.lstyle = {'FontSize', 8, 'BackgroundColor', c, 'EdgeColor', 0.8*c};
theme.plot_sol = @default_plot_sol;

if nargin>=1
  TB = tb_info.tb;
  ST = [tb_info.fam tb_info.(TB).branch_type];
  thmfunc = str2func(sprintf('%s_plot_theme', tb_info.tb));
  thm = thmfunc(ST);
  thm.special = union(theme.special, thm.special);
  theme = coco_merge(theme, thm);
end
end

function default_plot_sol(varargin)
error('%s: no plot function defined for plotting solutions\n');
end
