function theme = po_plot_theme(BT)
%PO_PLOT_THEME    Default plot theme for PO toolbox.
%
%   PO_PLOT_THEME with no arguments displays a list of default themes for
%   all branch types of the PO toolbox.
%
%   PO_PLOT_THEME(BT) returns a struct defining the plot theme for a branch
%   of type BT. This theme struct contains only settings for PO that
%   override or augment the default COCO plot theme. The full default theme
%   struct is computed in COCO_PLOT_THEME as
%   COCO_MERGE(COCO_PLOT_THEME,PO_PLOT_THEME(BT)).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: po_plot_theme.m 3063 2017-11-02 20:22:06Z hdankowicz $

if nargin<1
  display_themes();
else
  
  theme = struct();
  theme.sol.RO   = {'k-', 'LineWidth', 1, 'Marker', '.'};
  theme.plot_sol = @po_plot_sol;
  
  switch BT
    
    case 'po'
      theme.bd.col2  = '||x||_{2,MPD}';
      lspec_s = {'b-', 'LineWidth', 1};
      lspec_u = {'b--', 'Linewidth', 1};
      theme.lspec    = {lspec_s, lspec_u};
      theme.ustab    = 'po.test.USTAB';
      theme.ustabfun = @(x) (x>=1)+1;
      theme.usept    = {'SN', 'PD', 'TR', 'FP', 'BP'};
      theme.SN       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8};
      theme.PD       = {'kd', 'MarkerFaceColor', 'r', 'MarkerSize', 8};
      theme.TR       = {'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8};
      theme.NSA      = {'kd', 'MarkerFaceColor', 'm', 'MarkerSize', 8};
      theme.sol.SN   = {'g-', 'LineWidth', 2};
      theme.sol.PD   = {'r-', 'LineWidth', 2};
      theme.sol.TR   = {'b-', 'LineWidth', 2};
      
    case 'po.SN'
      theme.lspec   = {'g-', 'LineWidth', 2};
      
    case 'po.PD'
      theme.lspec   = {'m-', 'LineWidth', 2};
      
    case 'po.TR'
      theme.lspec   = {'c-', 'LineWidth', 2};
      
    otherwise
      error('%s: unknown solution branch type ''%s''', mfilename, ST);
  end

end

end

function display_themes()
tb_info.tb = 'po';

fprintf('PO default plotting themes:\n');

tb_info.po.branch_type = 'po';
default = coco_plot_theme(tb_info);
disp(' ');
disp('''po'' =');
disp(' ');
disp(default);

tb_info.po.branch_type = 'po.SN';
SN = cocoe_plot_theme(tb_info);
disp(' ');
disp('''po.SN'' =');
disp(' ');
disp(SN);

tb_info.po.branch_type = 'po.PD';
PD = coco_plot_theme(tb_info);
disp(' ');
disp('''po.PD'' =');
disp(' ');
disp(PD);

tb_info.po.branch_type = 'po.TR';
TR = coco_plot_theme(tb_info);
disp(' ');
disp('''po.TR'' =');
disp(' ');
disp(TR);

end
