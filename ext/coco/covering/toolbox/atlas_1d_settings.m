function atlas_1d_settings(varargin)
%ATLAS_1D_SETTINGS   Show and explain settings of 'atlas_1d' atlas algorithm.
%
% ATLAS_1D_SETTINGS(VARARGIN)
%
% VARARGIN = { }
% Show and explain default settings.
%
% VARARGIN = { PROB }
% Show active settings.
%
% PROB : Continuation problem structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: atlas_1d_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

str = coco_stream(varargin{:});
if isstruct(str.peek)
  prob = str.get;
else
  prob = coco_prob;
end

switch numel(str)
  
  case 0
    cont = coco_get(prob, 'cont');
    [sets,spec] = atlas_1d.get_settings(cont);
    help atlas_1d_settings>main_settings
    if nargin==0
      coco_explain_settings(spec, sets);
      fprintf('\n');
      atlas_settings(1);
    else
      coco_explain_settings(spec, sets);
      fprintf('\n');
      atlas_settings(prob,1);
    end
    
  otherwise
    help atlas_1d_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of 'atlas_1d' atlas algorithm:
%============================
%
%PtMX   : Maximum # of steps in either direction along the solution manifold
%NAdapt : # of steps between each adaptive remeshing
%h0     : Initial continuation step size
%h_max  : Maximum continuation step size
%h_min  : Minimum continuation step size
%
%Full list of all 'atlas_1d' settings with default values:
%-------------------------------------------------
end
