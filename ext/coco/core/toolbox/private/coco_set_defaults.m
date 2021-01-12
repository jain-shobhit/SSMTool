function [opts] = coco_set_defaults(varargin)
% COCO_SET_DEFAULTS  Initialise COCO options for continuation.
%
%   OPTS = COCO_SET_DEFAULTS()
%   OPTS = COCO_SET_DEFAULTS(OPTS) sets a number of default property values
%   of the options structure. Type 'help coco_opts' to see a listing of the
%   various classes and properties set by COCO_SET_DEFAULTS. This function
%   is invoked by COCO.
%
%   See also: coco, coco_set, coco_opts.
%

error('%s: function obsolete', mfilename);

if nargin<=0
	opts = [];
else
	opts = varargin{1};
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct default hirarchy of functions

opts = coco_set(opts, 'pdat', 'vectorised', 'true');

opts = coco_set(opts, 'efunc', 'F',          @efunc_F);
opts = coco_set(opts, 'efunc', 'DFDX',       @efunc_DFDX);
opts = coco_set(opts, 'efunc', 'monitor_F',  @efunc_monitor_F);
opts = coco_set(opts, 'efunc', 'events_F',   @efunc_events_F);
opts = coco_set(opts, 'efunc', 'init',       @efunc_init);
opts = coco_set(opts, 'efunc', 'linsolve',   @default_linsolve);
% opts = coco_set(opts, 'efunc', 'func',     ); % must be set by toolboxes

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set defaults for corrector (newton's method)

opts = coco_set(opts, 'nwtn', 'ItMX',     10    ); % max. number of iterations
opts = coco_set(opts, 'nwtn', 'ItNW',     []    ); % max. number of Newton iterations
opts = coco_set(opts, 'nwtn', 'SubItMX',  8     ); % max. number of damping steps
opts = coco_set(opts, 'nwtn', 'TOL',      1.0e-8); % convergence criterion
opts = coco_set(opts, 'nwtn', 'LogLevel', 1     ); % level of diagnostic output
opts = coco_set(opts, 'nwtn', 'MaxStep',  0.1   ); % max. size of Newton step
opts = coco_set(opts, 'nwtn', 'ga0',      1.0   ); % initial damping factor
opts = coco_set(opts, 'nwtn', 'al',       0.5   ); % increase damping factor

opts = coco_set(opts, 'nwtn', 'func',     'func'); % functions to use

opts = coco_set(opts, 'nwtn', 'print_headline', @default_print);
opts = coco_set(opts, 'nwtn', 'print_data',     @default_print);

opts = coco_set(opts, 'nwtn', 'begin_callback',      @default_callback);
opts = coco_set(opts, 'nwtn', 'begin_step_callback', @default_callback);
opts = coco_set(opts, 'nwtn', 'end_callback',        @default_callback);
opts = coco_set(opts, 'nwtn', 'end_step_callback',   @default_callback);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set default covering toolbox

if exist('covering_create', 'file')
  opts = coco_set(opts, 'all', 'ContAlg',   'covering');
else
  opts = coco_set(opts, 'all', 'ContAlg',   @cover1d);
end
opts = coco_set(opts, 'all', 'callback',  @default_callback);
opts = coco_set(opts, 'all', 'linsolve',  @default_linsolve);
opts = coco_set(opts, 'all', 'print',     @default_print);
opts = coco_set(opts, 'all', 'save',      @default_save);
opts = coco_set(opts, 'all', 'update',    @default_update);
opts = coco_set(opts, 'all', 'CleanData', 0);
