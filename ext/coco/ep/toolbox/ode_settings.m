function ode_settings(varargin)
%ODE_SETTINGS   Show and explain settings of toolbox family ODE.
%
% ODE_SETTINGS(VARARGIN)
%
% This function is part of the toolbox developers interface and should
% typically not be used directly.
%
% VARARGIN = { }
% Show and explain default settings of ODE.
%
% VARARGIN = { PROB }
% Show active settings of ODE.
%
% VARARGIN = { PROB OID }
% Show active settings of instance of ODE appended to OID in PROB.
%
% VARARGIN = { PROB OID DATA }
% Show settings of ODE that will be active when restarting from a saved
% solution with ODE toolbox family data structure DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ode_settings.m 3157 2019-12-16 23:36:58Z hdankowicz $

str = coco_stream(varargin{:});
if isstruct(str.peek)
  prob = str.get;
else
  prob = coco_prob;
end

switch numel(str)
  
  case 0
    ode = struct;
    tbid = 'ode';
%     [~,spec] = ode_get_settings(prob, tbid, ode);
    [sets,spec] = ode_get_settings(prob, tbid, ode);
    if nargin==0
      help ode_settings>main_settings
    end
%     coco_explain_settings(prob, spec, ode, tbid);
    coco_explain_settings(spec, sets);
    
  case 1
    ode = struct;
    tbid = coco_get_id(str.get, 'ode');
    %     [~,spec] = ode_get_settings(prob, tbid, ode);
    %     coco_explain_settings(prob, spec, ode, tbid);
    [sets,spec] = ode_get_settings(prob, tbid, ode);
    coco_explain_settings(spec, sets);
    
  case 2
    tbid = coco_get_id(str.get, 'ode');
    data = str.get;
%     [~,spec] = ode_get_settings(prob, tbid, data.ode);
%     coco_explain_settings(prob, spec, data.ode, tbid);
    [sets,spec] = ode_get_settings(prob, tbid, data.ode);
    coco_explain_settings(spec, sets);
    
  otherwise
    help ode_settings
    
end
end

function main_settings %#ok<DEFNU>
%Main settings of ODE toolbox:
%=============================
%
%vectorized : The vectorized property enables a reduction in the number of
%  function evaluations required to compute all the columns of the Jacobian
%  matrix of f, and might significantly reduce computation time. A
%  vectorized evaluation of f(x,p) is possible if f is encoded such that
%  f([x1 x2 ...], [p1 p2 ...]) returns [f(x1,p1) f(x2,p2), ...]. This
%  enables fast evaluation of finite differences.
%
%  Set to off/false to disable vectorized evaluation. See ODESET for more
%  details on vectorized evaluation of functions.
%
%autonomous : The autonomous property indicates whether a given
%  ODE is autonomous or not. An autonomous ode has a time-independent
%  right-hand side f as in
%
%    x' = dx/dt = f(x,p),
%
%  whereas a non-autonomous ODE has a time-dependent right-hand side f as
%  in
%
%    x' = dx/dt = f(t,x,p).
%
%  Since dynamical systems are commonly defined by autonomous ODEs, this
%  setting is true by default. However, a number of toolboxes support
%  continuation of solutions to non-autonomous ODEs as well. Set autonomous
%  to off/false to change to time-dependent evaluation for toolboxes that
%  support non-autonomous problems.
%
%Full list of all ODE settings with default values:
%--------------------------------------------------
end
