function theme = hspo_plot_theme(BT)
%HSPO_PLOT_THEME    Default plot theme for HSPO toolbox.
%
%   HSPO_PLOT_THEME with no arguments displays a list of default themes for
%   all branch types.
%
%   HSPO_PLOT_THEME(BT) returns a struct defining the plot theme for a
%   branch of type BT. This theme struct contains only settings for HSPO
%   that override or augment the default COCO plot theme. The full default
%   theme struct is computed in COCO_PLOT_THEME as
%   COCO_MERGE(COCO_PLOT_THEME,HSPO_PLOT_THEME(BT)).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_plot_theme.m 2948 2016-04-10 10:41:39Z fschild $

if nargin<1
  display_themes();
else
  
  theme = struct();
  theme.sol.RO   = {'k-', 'LineWidth', 1, 'Marker', '.'};
  theme.plot_sol = @hspo_plot_sol;
  
  switch BT
    
    case 'hspo'
      lspec_s = {'b-', 'LineWidth', 1};
      lspec_u = {'b--', 'Linewidth', 1};
      theme.lspec    = {lspec_s, lspec_u};
      theme.ustab    = 'hspo.test.USTAB';
      theme.ustabfun = @(x) (x>=1)+1;
      theme.usept    = {'SN', 'PD', 'TR', 'FP', 'BP'};
      theme.SN       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8};
      theme.PD       = {'kd', 'MarkerFaceColor', 'r', 'MarkerSize', 8};
      theme.TR       = {'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8};
      theme.NS       = {'kd', 'MarkerFaceColor', 'm', 'MarkerSize', 8};      
      theme.sol.SN   = {'g-', 'LineWidth', 2};
      theme.sol.PD   = {'r-', 'LineWidth', 2};
      theme.sol.TR   = {'b-', 'LineWidth', 2};
      
    case 'hspo.SN'
      theme.lspec   = {'g-', 'LineWidth', 2};
      
    case 'hspo.PD'
      theme.lspec   = {'m-', 'LineWidth', 2};
      
    case 'hspo.TR'
      theme.lspec   = {'c-', 'LineWidth', 2};
      
    otherwise
      error('%s: unknown solution branch type ''%s''', mfilename, BT);
  end

end

end

function display_themes()
tb_info.tb = 'hspo';

fprintf('HSPO default plotting themes:\n');

tb_info.po.branch_type = 'hspo';
default = ode_plot_theme(tb_info);
disp(' ');
disp('''hspo'' =');
disp(' ');
disp(default);

tb_info.po.branch_type = 'hspo.SN';
SN = ode_plot_theme(tb_info);
disp(' ');
disp('''hspo.SN'' =');
disp(' ');
disp(SN);

tb_info.po.branch_type = 'hspo.PD';
PD = ode_plot_theme(tb_info);
disp(' ');
disp('''hspo.PD'' =');
disp(' ');
disp(PD);

tb_info.po.branch_type = 'hspo.TR';
NS = ode_plot_theme(tb_info);
disp(' ');
disp('''hspo.TR'' =');
disp(' ');
disp(TR);

end
