function [cont, spec] = get_settings(cont)

spec = {% Name    Type   Default  Action  Args  Description
        'PtMX',  '[int]',     100,   'read', {}, 'max # of continuation steps'
      'NAdapt',    'int',       0,   'read', {}, 'Adaptation period, 0 == off'
          'h0',    'num',     0.1,   'read', {}, 'initial step size'
       'h_max',    'num',     0.5,   'read', {}, 'max step size'
       'h_min',    'num',    0.01,   'read', {}, 'min step size'
          'FP', 'log|on',    true, 'switch', {}, 'detect fold points'
        'fpar',    'str',      '',   'read', {}, 'active continuation parameter for fold detection'
          'BP', 'log|on',    true, 'switch', {}, 'detect branch points'
        'RMMX',    'int',      10,   'read', {}, 'max # of remesh sweeps'
   'h_fac_max',    'num',       2,   'read', {}, 'max step size adaptation factor'
   'h_fac_min',    'num',     0.5,   'read', {}, 'min step size adaptation factor'
      'MaxRes',    'num',     0.1,   'read', {}, 'max residual norm in prediction'
      'al_max',    'num',       7,   'read', {}, 'max angle between consecutive tangents'
          'ga',    'num',    0.95,   'read', {}, 'adaptation security factor'
   'bi_direct', 'log|on',    true, 'switch', {}, 'go in both directions or not'
      'interp',    'str', 'cubic',   'read', {}, 'cseg interpolation'
      'Valpha',    'num',      80,   'read', {}, 'tolerance for "vertical" tangent'
    'NullItMX',    'int',       0,   'read', {}, 'max # of nullspace corrections'
  };

cont = coco_parse_settings([], spec, cont, '');

cont.h0        = min(max(cont.h0, cont.h_min), cont.h_max);
cont.arc_alpha = cont.al_max * pi / 180;

if isfield(cont, 'ItMX')
  cont.PtMX = cont.ItMX;
end
assert(any(numel(cont.PtMX)==[1 2]));
if numel(cont.PtMX)==1
  if cont.bi_direct
    cont.PtMX = [ cont.PtMX cont.PtMX ];
  else
    if cont.PtMX>=0
      cont.PtMX = [ 0 cont.PtMX ];
    else
      cont.PtMX = [ cont.PtMX 0 ];
    end
  end
end
cont.PtMX = [-1 1] .* abs(cont.PtMX);

end
