function install
% run this file first to install all external packages  and
% also the main software appropriately

maindir = fileparts(mfilename('fullpath'));
extdir = 'ext';
cocoinstall = fullfile(maindir, extdir , 'coco','startup.m');
run(cocoinstall);

NLvibpath = fullfile(maindir, extdir, 'NLvib_v1.3','SRC');
addpath(genpath(NLvibpath))

addpath(fullfile(maindir, extdir,'combinator'));

addpath(fullfile(maindir, extdir, 'tensor_toolbox'));

addpath(fullfile(maindir, extdir, 'Wrappers'));

addpath(fullfile(maindir, 'src'));

addpath(fullfile(maindir, 'src','misc'));

addpath(fullfile(maindir, 'src','multiindex'));

addpath(fullfile(maindir, 'src', 'frc'));
end


