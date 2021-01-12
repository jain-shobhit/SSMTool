function [opts sol] = save_data(opts, pt1, pt0)
%SAVE_DATA   Save point data.
%
%   OPTS = COCO_SAVE_DATA(OPTS) is called whenever a new solution has been
%   computed. COCO_SAVE_DATA checks if a label was or should be assigned to
%   the current solution and calls the appropriate save function if so. It
%   inserts the current solution point into the bifurcation diagram.
%
%   See also: save_full, save_reduced, bddat_insert_point
%

%% we might modify opts, so check that opts gets assigned

if nargout<1
	error('%s: too few output arguments', mfilename);
end

%% save solution point

sol      = pt1;
sol.lab  = [];
next_lab = opts.bddat.next_lab;
NSV      = opts.cont.NSV;
if isempty(NSV)
	NSV    = opts.cont.NPR;
end

if nargin >= 3
	
	% MX point (correct failed), full save
  if isempty(pt0.pt_type)
    sol.pt_type = 'MX';
  else
    sol.pt_type = pt0.pt_type;
  end
	sol.lab     = next_lab;
	next_lab    = next_lab + 1;
	save_func   = @coco_save_full;

elseif isfield(sol, 'pt_type') && ~isempty(sol.pt_type)
	
	% special solution point, full save
	sol.lab   = next_lab;
	next_lab  = next_lab + 1;
	save_func = @coco_save_full;

elseif mod(sol.pt, NSV) == 0
	
	% regular output point for restart, full save
	sol.pt_type = 'RO';
	sol.lab     = next_lab;
	next_lab    = next_lab + 1;
	save_func   = @coco_save_full;
	
elseif mod(sol.pt, opts.cont.NPR) == 0
	
	% regular output point for plotting, reduced save
	sol.pt_type = 'ROS';
	sol.lab     = next_lab;
	next_lab    = next_lab + 1;
	save_func   = @coco_save_reduced;
	
end

opts.bddat.next_lab = next_lab;

%% insert point into opts.bd and save bifurcation diagram

opts = opts.bddat.insert(opts, sol);

if ~isempty(sol.lab)
	if nargin>2
		opts = save_func(opts, sol, pt0);
	else
		opts = save_func(opts, sol);
  end
  opts = coco_bd_save(opts);
end
