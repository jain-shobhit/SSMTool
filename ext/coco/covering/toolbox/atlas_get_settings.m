function [cont, spec] = atlas_get_settings(cont, dim)

spec = {% Name    Type   Default  Action  Args  Description
        'NPR', 'int',      10,   'read', {}, 'frequency of screen outputs'
        'NSV', 'int',      10,   'read', {}, 'frequency of storing solutions to disk (default to ''NPR'')'
  'corrector', 'str',  'nwtn',   'read', {}, 'nonlinear corrector'
   'linsolve', 'str',  'splu',   'read', {}, 'linear solver'
      'atlas', 'str',      '',   'read', {}, 'atlas algorithm suffix'
  };

cont = coco_parse_settings([], spec, cont, '');

if nargin==2 && isempty(cont.atlas)
  switch dim
    case 0
      cont.atlas = '0d';
    case 1
      cont.atlas = '1d';
    otherwise
      cont.atlas = 'kd';
  end
end

end
