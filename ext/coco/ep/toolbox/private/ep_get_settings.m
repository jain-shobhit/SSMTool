function [ep, spec] = ep_get_settings(prob, tbid, ep)
%EP_GET_SETTINGS   Merge EP toolbox settings with default values.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_get_settings.m 2861 2015-07-26 22:08:48Z hdankowicz $

spec = {% Name    Type   Default  Action  Args  Description
       'bifus', 'log|on',  true, 'switch', {}, 'enable/disable detection of bifurcations'
       'USTAB', 'log|on',  true, 'switch', {}, 'monitor number of unstable eigenvalues'
          'SN', 'log|on',  true, 'switch', {}, 'detect saddle-node bifurcations'
          'HB', 'log|on',  true, 'switch', {}, 'detect Hopf bifurcations'
         'NSA', 'log|on', false, 'switch', {}, 'detect neutral saddle points'
         'BTP', 'log|on',  true, 'switch', {}, 'detect Bogdanov-Takens points'
  };

ep = coco_parse_settings(prob, spec, ep, tbid);

end
