% Toolbox: 'atlas_kd'
%
% This covering toolbox implements a k-dimensional expanding boundary atlas
% algorithm with adaptive step size and theta-based projection condition
% and predictor that recognizes computational domain boundaries and may
% start from a boundary. The atlas stores a chart network to prevent
% redundant coverage, and maintains a record of a polygonal boundary to
% prevent premature termination. The algorithm adapts the step size to
% local properties of the manifold and allows for local remeshing
% associated with a change in radius of a previously merged chart.
%
% Class settings
%   R     - Initial step size (default: 0.1)
%   R_max - Maximum step size (default: 1.0)
%   R_min - Minimum step size (default: 0.001)
%   R_fac_min - Minimum step size adaptation factor (default: 0.5)
%   R_fac_max - Maximum step size adaptation factor (default: 2.0)
%   ga    - Adaptation security factor (default: 0.95)
%   MaxRes - Maximal residuum for prediction step (default: 0.1)
%   PtMX  - Maximum number of continuation steps (default: 50)
%   theta - Theta method (default: 0.5)
%   almax - Critical angle between successive tangent vectors (default: 35)
%   Rmarg - Margin for merging charts into boundary (default: 1.05)

% This file is part of the atlas_kd toolbox, copyright (C) Michael
% Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
% COCO (http://sourceforge.net/projects/cocotools).
