function coco_recipes_doc(varargin)
% Show tutorial documentation in help browser.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_recipes_doc.m 2839 2015-03-05 17:09:01Z fschild $

p = fileparts(mfilename('fullpath'));
p = fullfile(p, varargin{:});
while ~isdir(p)
  p = fileparts(p);
end
p = addpath(p);
doc(varargin{end});
path(p);
end
