function prob = lsol_cond(prob)
%LSOL_COND   Append the condition number test function
%
% The constructor can only be called after all embedded monitor functions
% have been added to the continuation problem structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: lsol_cond.m 2955 2017-01-10 15:19:26Z hdankowicz $

prob = coco_add_func(prob, 'cond', @lsol_cond_TF, [], ...
  'regular', 'lsol.cond', 'PassChart', 'fdim', 1);
end

function [data chart y] = lsol_cond_TF(prob, data, chart, u)

cdata = coco_get_chart_data(chart, 'lsol'); % Extract chart data
if isfield(cdata, 'cond') % Linear solver computes this field
  y = cdata.cond;
else
  y = nan;
end

end
