function atlas_settings(varargin)
%ATLAS_SETTINGS   Show and explain settings of coco-compatible atlas algorithms.
%
% ATLAS_SETTINGS(VARARGIN)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% VARARGIN = { }
% Show and explain default settings of coco-compatible atlas algorithms.
%
% VARARGIN = { PROB }
% Show active settings of coco-compatible atlas algorithms.
%
% VARARGIN = { DIM }
% Show and explain default settings of coco-compatible atlas algorithms of dimension DIM.
%
% VARARGIN = { PROB DIM }
% Show and explain default settings of coco-compatible atlas algorithms of dimension DIM.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

str = coco_stream(varargin{:});
if isstruct(str.peek)
  prob = str.get;
else
  prob = coco_prob;
end

switch numel(str)
  
  case 0
    cont = coco_get(prob, 'cont');
    [sets,spec] = atlas_get_settings(cont);
    help atlas_settings>main_settings
    coco_explain_settings(spec, sets);
    
  case 1
    dim = str.get;
    cont = coco_get(prob, 'cont');
    [sets,spec] = atlas_get_settings(cont,dim);
    help atlas_settings>main_settings
    coco_explain_settings(spec, sets);

  otherwise
    help atlas_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of COCO-compatible atlas algorithms:
%=============================
%
%NPR   : Frequency of screen outputs
%NSV   : Frequency of storing solutions to disk (defaults to NPR)
%atlas : atlas algorithm suffice (e.g., '1d')
%
%Full list of all COCO-compatible atlas algorithm settings with default values:
%--------------------------------------------------
end
