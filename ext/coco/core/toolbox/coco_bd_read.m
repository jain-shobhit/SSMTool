function varargout = coco_bd_read(runid, varargin)

bdfname = coco_fname(runid, 'bd.mat');

if numel(varargin)==0
  list = { 'bd' };
else
  list = varargin;
end

vars = who('-file', bdfname);

if any(strcmp('version', vars))
  load(bdfname, 'version');
  switch version
    case 1
      [varargout{1:nargout}] = read_bd_v1(bdfname, list);
    otherwise
      error('%s: unknown format of bd file', mfilename);
  end
else
  if any(strcmp('bd_data', vars))
    [varargout{1:nargout}] = read_bd_v1(bdfname, list);
  else
    varargout{1} = read_bd_v0(bdfname);
  end
end

end

function varargout = read_bd_v1(bdfname, list)
S = load(bdfname, 'bd_data');
for i=1:max(1,nargout)
  varargout{i} = coco_slot_data(list{i}, S.bd_data);
end
end

function bd = read_bd_v0(bdfname)
S = load(bdfname, 'bd');
bd = S.bd;
end
