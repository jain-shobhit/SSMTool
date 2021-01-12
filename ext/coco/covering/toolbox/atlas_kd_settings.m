function atlas_kd_settings(varargin)
%ATLAS_KD_SETTINGS   Show and explain settings of 'atlas_kd' atlas algorithm.
%
% ATLAS_KD_SETTINGS(VARARGIN)
%
% VARARGIN = { DIM }
% Show and explain default settings.
%
% VARARGIN = { PROB DIM }
% Show active settings.
%
% PROB : Continuation problem structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: atlas_1d_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

str = coco_stream(varargin{:});
if isstruct(str.peek)
  prob = str.get;
  assert(nargin==2, '%s: missing manifold dimension', mfilename);
  dim  = str.get;
else
  prob = coco_prob;
  assert(nargin==1, '%s: missing manifold dimension', mfilename);
  dim  = str.get;
end

switch numel(str)
  
  case 0
    cont = coco_get(prob, 'cont');
    [sets,spec] = atlas_kd.get_settings(cont);
    help atlas_kd_settings>main_settings
    if nargin==0
      coco_explain_settings(spec, sets);
      fprintf('\n');
      atlas_settings(dim);
    else
      coco_explain_settings(spec, sets);
      fprintf('\n');
      atlas_settings(prob,dim);
    end
    
  otherwise
    help atlas_kd_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of 'atlas_kd' atlas algorithm:
%============================
%
%PtMX   : Maximum # of steps along the solution manifold
%NAdapt : # of steps between each adaptive remeshing
%R      : Initial continuation step size
%R_max  : Maximum continuation step size
%R_min  : Minimum continuation step size
%
%Full list of all 'atlas_kd' settings with default values:
%-------------------------------------------------
end
