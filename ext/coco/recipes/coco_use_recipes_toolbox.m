function coco_use_recipes_toolbox( varargin )
%COCO_USE_RECIPES_TOOLBOX   Add or remove tutorial toolboxes to search path.
%
% varargin = sequence of toolbox folder names
%
% When nargin = 0, function executes the coco_rm_this_path function
% in every folder on the search path.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_use_recipes_toolbox.m 2839 2015-03-05 17:09:01Z fschild $

p   = fileparts(mfilename('fullpath'));
tbx = {};

for i = 1:nargin
  t = fullfile(p, varargin{i});
  assert(exist(t, 'dir')==7, 'toolbox ''%s'' not found', t);
  tbx = [ tbx t ]; %#ok<AGROW>
end

while exist('coco_rm_this_path', 'file')==2
  coco_rm_this_path;
end

if nargin>0
  addpath(tbx{:});
end

end
