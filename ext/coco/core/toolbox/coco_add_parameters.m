function opts = coco_add_parameters(varargin)
%COCO_ADD_PARAMETERS   Add external parameters to continuation problem.
%
% Obsolete. See also: COCO_ADD_PARS.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_parameters.m 2839 2015-03-05 17:09:01Z fschild $

fprintf(2, '%s: function will become obsolete, use ''coco_add_pars'' instead.\n', ...
  mfilename);

opts = coco_add_pars(varargin{:});

end

