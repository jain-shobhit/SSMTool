function coco_recipes_edit(varargin)
% Open tutorial file in editor.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_recipes_edit.m 2839 2015-03-05 17:09:01Z fschild $

p = fileparts(mfilename('fullpath'));
p = fullfile(p, varargin{:});

[mfile remain] = strtok(p, '>');
if isempty(remain)
 edit(p);
else
 func = strtok(remain, '>');
 % bug: replace use of unsupported functions
 if ~isempty(which('editorservices.openAndGoToFunction'))
   editorservices.openAndGoToFunction([mfile '.m'], func);
 elseif ~isempty(which('matlab.desktop.editor.openAndGoToFunction'))
   matlab.desktop.editor.openAndGoToFunction([mfile '.m'], func);
 else
   edit(mfile); % fall back
 end
end
