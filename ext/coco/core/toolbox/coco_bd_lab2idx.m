function idx = coco_bd_lab2idx(bd, labs)
%COCO_BD_LAB2IDX    Compute indices corresponding to labels.
%
%   Usage:
%
%   idx = coco_bd_lab2idx(bd, 1);   % get row-index of label 1
%   x   = coco_bd_col(bd, '||U||'); % extract column ||U||
%   y   = x(idx);                   % get ||U|| for solution with label 1

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_bd_lab2idx.m 2946 2016-04-07 12:54:52Z fschild $

LCOL = coco_bd_col(bd, 'LAB', false);
idx  = zeros(size(labs));
N    = numel(labs);

for i=1:N
  func = @(x) ~isempty(x) && (x==labs(i));
  idx(i) = find(cellfun(func, LCOL), 1, 'first');
end
end
