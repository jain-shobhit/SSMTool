function [coll, spec] = coll_get_settings(prob, tbid, coll)
%COLL_GET_SETTINGS   Merge COLL toolbox settings with default values.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_get_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

if ~coco_exist('TOL', 'class_prop', prob, tbid, '-no-inherit-all')
  TOL = coco_get(prob, 'corr', 'TOL')^(2/3);
  prob = coco_set(prob, tbid, 'TOL', TOL);
end
TOL = coco_get(prob, 'coll', 'TOL');

spec1 = {% Name   Type   Default   Action   Args   Description
         'var', 'log|on', false, 'switch', {}, 'enable/disable temporary storage of solution to variational problem'
        'NTST',    'int',    10,   'read', {}, 'number of mesh intervals'
        'NCOL',    'int',     4,   'read', {}, 'number of collocation nodes'
         'SAD',    'num',  0.95,   'read', {}, 'equidistribution weight for error estimator'
      'method',    'str',  '3I',   'read', {}, 'choice of Banach iteration boundary condition'
       'NBeta',    'int',     5,   'read', {}, 'number of homotopy steps for initialization of fundamental solution'
      'NBItMX',    'int',    10,   'read', {}, 'maximum number of Banach iterations for fundamental solution'
         'TOL',    'num',   TOL,   'read', {}, 'discretization error tolerance'
        'MXCL', 'log|on',  true, 'switch', {}, 'enable/disable termination when discretization error exceeds tolerance'
  };
coll = coco_parse_settings(prob, spec1, coll, tbid);

assert(coll.SAD<=1 && coll.SAD>=0, ...
  '%s: input for option ''SAD'' not in [0,1]', tbid);

spec2 = {% Name   Type   Default   Action   Args   Description
   'TOLINC', 'num',          coll.TOL/5, 'read', {}, 'Upper bound on discretization error in window of adaptation'
   'TOLDEC', 'num',         coll.TOL/20, 'read', {}, 'Lower bound on discretization error in window of adaptation'
   'NTSTMN', 'int', min(  5, coll.NTST), 'read', {}, 'Minimum number of discretization intervals'
   'NTSTMX', 'int', max(100, coll.NTST), 'read', {}, 'Maximum number of discretization intervals'
  };
coll = coco_parse_settings(prob, spec2, coll, tbid);

spec = [spec1; spec2];

end
