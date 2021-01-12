function cp = coco_path()
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
cp   = coco_genpaths(root);
end

function cp = coco_genpaths(root) % bug: update to new repo layout.
tbroot = fullfile(root, 'toolbox');
exroot = fullfile(root, 'examples');

cp = {};

exas = { fullfile(exroot, 'recipes', 'example') };
for i=1:numel(exas)
  exdir = exas{i};
  if exist(exdir, 'dir')
    cp = [ cp ; { exdir } ]; %#ok<AGROW>
  end
end

tbxs = dir(tbroot);
if isempty(tbxs)
  return
end
isdir = logical(cat(1,tbxs.isdir));
tbxs  = { tbxs(isdir).name };
tbxs  = tbxs(~strncmp('.', tbxs, 1));
for i=1:numel(tbxs)
  tbdir = fullfile(tbroot, tbxs{i}, 'toolbox');
  cp    = [cp ; coco_genpath(tbdir, {'private'})]; %#ok<AGROW>
end

end

function p = coco_genpath(d, ignore)

p = {d};

files = dir(d);
if isempty(files)
  return
end

isdir = logical(cat(1,files.isdir));
dirs  = files(isdir);
prune = { '.' '@' '+' };

for i=1:length(dirs)
  dirname = dirs(i).name;
  if ~( any(strncmp(dirname, prune, 1)) || any(strcmp(dirname, ignore)) )
    p = [p ; coco_genpath(fullfile(d,dirname), ignore)]; %#ok<AGROW>
  end
end

end
