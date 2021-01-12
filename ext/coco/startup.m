function startup
cocodir = fileparts(mfilename('fullpath'));
addpath(fullfile(cocodir, 'core', 'toolbox'));
addpath(fullfile(cocodir, 'covering', 'toolbox'));
addpath(fullfile(cocodir, 'ep', 'toolbox'));
addpath(fullfile(cocodir, 'coll', 'toolbox'));
addpath(fullfile(cocodir, 'po', 'toolbox'));
addpath(fullfile(cocodir, 'recipes'));
addpath(fullfile(cocodir, 'continex', 'toolbox'));
end

% function startup
% cocodir = fileparts(fileparts(mfilename('fullpath')));
% addpath(fullfile(cocodir, 'toolbox'));
% cp = coco_path;
% addpath(cp{:});
% end
