function [hspo, spec] = hspo_get_settings(prob, tbid, hspo)
%PHSO_GET_SETTINGS   Merge HSPO toolbox settings with default values.

spec = {% Name    Type   Default  Action  Args  Description
       'bifus', 'log|on',  true, 'switch', {}, 'enable/disable detection of bifurcations'
       'USTAB', 'log|on',  true, 'switch', {}, 'monitor number of unstable eigenvalues'
          'SN', 'log|on',  true, 'switch', {}, 'detect saddle-node bifurcations'
          'PD', 'log|on',  true, 'switch', {}, 'detect period-doubling bifurcations'
          'TR', 'log|on',  true, 'switch', {}, 'detect torus bifurcations'
         'NSA', 'log|on', false, 'switch', {}, 'detect neutral saddle points'
  };

hspo = coco_parse_settings(prob, spec, hspo, tbid);

end
