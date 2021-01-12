function theme = coll_plot_theme(BT)
%COLL_PLOT_THEME    Default plot theme for COLL toolbox.
%
%   COLL_PLOT_THEME with no arguments displays a list of default themes for
%   all branch types of the COLL toolbox.
%
%   COLL_PLOT_THEME(BT) returns a struct defining the plot theme for a
%   branch of type BT. This theme struct contains only settings for COLL
%   that override or augment the default COCO plot theme. The full default
%   theme struct is computed in COCO_PLOT_THEME as
%   COCO_MERGE(COCO_PLOT_THEME, COLL_PLOT_THEME(BT)).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_plot_theme.m 2948 2016-04-10 10:41:39Z fschild $

if nargin<1
  display_themes();
else
  
  theme = struct();
  theme.bd.col2  = '||x||_{L_2[0,1]}';
  theme.lspec    = {'k-', 'LineWidth', 1};
  theme.sol.col1 = 't';
  theme.sol.col2 = 'x';
  theme.sol.RO   = {'k-', 'LineWidth', 2, 'Marker', '.', 'MarkerSize', 12};
  theme.plot_sol = @coll_plot_sol;
  
  switch BT
    
    case {'seg', 'seg.VAR'}
      
    otherwise
      error('%s: unknown solution branch type ''%s''', mfilename, BT);
  end
  
end

end

function display_themes()
tb_info.tb = 'coll';

fprintf('COLL default plotting themes:\n');

tb_info.coll.branch_type = 'seg';
default = ode_plot_theme(tb_info);
disp(' ');
disp('''seg'' =');
disp(' ');
disp(default);

tb_info.coll.branch_type = 'seg.var';
default = ode_plot_theme(tb_info);
disp(' ');
disp('''seg.var'' =');
disp(' ');
disp(default);

end
