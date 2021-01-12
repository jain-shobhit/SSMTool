function [po, spec] = po_get_settings(prob, tbid, po)
%PO_GET_SETTINGS   Merge PO toolbox settings with default values.

spec = {% Name    Type   Default  Action  Args  Description
       'bifus', 'log|on',  true, 'switch', {}, 'enable/disable detection of bifurcations'
       'USTAB', 'log|on',  true, 'switch', {}, 'monitor number of unstable eigenvalues'
          'SN', 'log|on',  true, 'switch', {}, 'detect saddle-node bifurcations'
          'PD', 'log|on',  true, 'switch', {}, 'detect period-doubling bifurcations'
          'TR', 'log|on',  true, 'switch', {}, 'detect torus bifurcations'
         'NSA', 'log|on', false, 'switch', {}, 'detect neutral saddle points'
  };

po = coco_parse_settings(prob, spec, po, tbid);

end
