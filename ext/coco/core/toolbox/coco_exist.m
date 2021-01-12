function flag = coco_exist(name, type, varargin)
% COCO_EXIST Check if coco-objects are defined.
%
% FLAG = COCO_EXIST(NAME, 'run')
%
% FLAG = COCO_EXIST(NAME, 'class_prop', PROB, PATH, [INHERIT-MODE])
% FLAG = COCO_EXIST(NAME, 'opts_path' , PROB, [INHERIT-MODE])
% FLAG = COCO_EXIST(NAME, 'func'      , PROB)
% FLAG = COCO_EXIST(NAME, 'slot'      , PROB)
% FLAG = COCO_EXIST(NAME, 'signal'    , PROB)
%
% FLAG = COCO_EXIST(NAME, 'col', BD)

switch lower(type)
  
  case 'run'
    % varargin = { }
    [fname flag] = coco_fname(name, 'bd.mat'); %#ok<ASGLU>
    
  case { 'class_prop' 'class_property' }
    % varargin = { prob path [inherit-mode] }
    opts = varargin{1}.opts;
    path = varargin{2};
    if nargin<5
      flag = opts.prop_exist(path, name, '-inherit');
    else
      flag = opts.prop_exist(path, name, varargin{3});
    end

  case 'opts_path'
    % varargin = { prob [inherit-mode] }
    opts = varargin{1}.opts;
    if nargin<4
      flag = opts.path_exist(name, '-inherit');
    else
      flag = opts.path_exist(name, varargin{2});
    end

  case { 'func' 'function' }
    % varargin = { prob }
    prob = varargin{1};
    flag = ( isfield(prob, 'efunc') && ...
      isfield(prob.efunc, 'identifyers') && ...
      any(strcmpi(name, prob.efunc.identifyers)) );
    
  case 'slot'
    % varargin = { prob }
    prob = varargin{1};
    flag = ( isfield(prob, 'slots') && isfield(prob.slots, lower(name)) );
    
  case { 'sig' 'signal' }
    % varargin = { prob }
    prob = varargin{1};
    flag = ( isfield(prob, 'signals') && isfield(prob.signals, lower(name)) );
    
  case { 'col' 'column' }
    % varargin = { bd }
    bd   = varargin{1};
    flag = ( any(strcmp(name, bd(1,:))) );
    
  otherwise
    if ischar(type)
      error('%s: unknown object type ''%s''', mfilename, type);
    else
      error('%s: unknown object type', mfilename);
    end
end

end
