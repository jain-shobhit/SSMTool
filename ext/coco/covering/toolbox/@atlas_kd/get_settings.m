function [cont, spec] = get_settings(cont, dim)
%GET_SETTINGS   Defines default class settings.
%
% Merge default class settings with cont argument.
%
% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).

spec = {% Name    Type   Default  Action  Args  Description
        'PtMX',  '[int]',     100,   'read', {}, 'max # of continuation steps'
      'NAdapt',    'int',       0,   'read', {}, 'Adaptation period, 0 == off'
        'RMMX',    'int',      10,   'read', {}, 'max # of remesh sweeps'   
           'R',    'num',     0.1,   'read', {}, 'initial step size'
       'R_max',    'num',       1,   'read', {}, 'max step size'
       'R_min',    'num',    0.01,   'read', {}, 'min step size'
   'R_fac_max',    'num',       2,   'read', {}, 'max step size adaptation factor'
   'R_fac_min',    'num',     0.5,   'read', {}, 'min step size adaptation factor'
          'ga',    'num',    0.95,   'read', {}, 'adaptation security factor'
      'MaxRes',    'num',     0.1,   'read', {}, 'max residual norm in prediction'
       'almax',    'num',      10,   'read', {}, 'max angle between consecutive tangents'
       'Rmarg',    'num',     1.1,   'read', {}, 'scaling factor for circumscribed cube'
       'theta',    'num',     0.5,   'read', {}, 'theta method parameter'
  };

cont = coco_parse_settings([], spec, cont, '');

cont.almax = cont.almax*pi/180;
if nargin==2
  cont.cube = create_cube(dim, cont.Rmarg*cont.R); % Initial dim-dimensional cube
  cont.s0   = cell2mat(cont.cube.v)/(cont.Rmarg*cont.R*sqrt(dim)); % Vertex directions
  cont.bv0  = 1:cont.cube.n; % Indices of available vertices
end

end
