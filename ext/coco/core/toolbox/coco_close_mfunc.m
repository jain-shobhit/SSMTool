function opts = coco_close_mfunc(opts)
%COCO_SET_DIM   Set dim. of manifold and define prim. cont. parameters.
%
%   [OPTS ARGNUM XP0] = COCO_SET_DIM(OPTS, [DIM,] CONT_PARS, PAR_WINDOW)
%   prepares the system EFUNC, composed of the zero problem and the monitor
%   functions, for continuation. After a successful call to COCO_SET_DIM
%   the function opts.efunc.F maps (N+DIM)-vectors to N-vectors, that is,
%   opts.efunc.F has a solution manifold of dimension DIM. The first DIM
%   parameters in CONT_PARS become the primary continuation parameters and
%   the remaining parameters are exchanged for internal parameters or used
%   for output during continuation. A typical manifold covering algorithm
%   will augment the system in EFUNC with DIM projection conditions to
%   obtain a unique solution point on the manifold.
%
%   OPTS is coco's options structure, DIM is the dimension of the solution
%   manifold and CONT_PARS is a cell array with at least DIM elements
%   containing the names of the continuation parameters.
%
%   See also: COCO_ADD_FUNC, COCO_SET_FUNC, COCO_ADD_EVENT

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coco_close_mfunc.m 3101 2019-06-07 17:33:01Z hdankowicz $

%% affected fields in opts
%
%    opts.efunc.dim           - dimension of solution manifold
%                              (copy of argument DIM)
%    opts.efunc.cont_pars     - list of continuation parameters
%                              (copy of argument CONT_PARS)
%    opts.efunc.cont_par_idx  - indices of continuation parameters
%    opts.efunc.userA_dim     - dimension of solution manifold of monitor
%                              functions (set to 0 if non-existent)
%
%    opts.efunc.acp_idx       - indices of active continuation parameters
%    opts.efunc.op_idx        - indices of output parameters
%
%    opts.efunc.ev.pidx       - indices of event parmeters
%    opts.efunc.ev.par_type   - type of event parameter
%    opts.efunc.ev.vals       - event values
%    opts.efunc.ev.point_type - event name (point type)
%    opts.efunc.ev.BP_idx     - indices of boundary events
%    opts.efunc.ev.MX_idx     - indices of terminate events
%    opts.efunc.ev.SP_idx     - indices of special point events

%% algorithm
%
%  1. create temporary copy of userA parameter indices in tmp
%  2. exchange parameters in tmp as specified in opts.efunc.xchg
%  3. compute indices of active continuation and output parameters
%  4. create event structure

%% add pending functions

opts = coco_add_pending(opts, 'mfunc');
if ~isempty(opts.efunc.pending)
  % bug: generate sensible error message
  flist = opts.efunc.pending{1,1};
  for i=2:size(opts.efunc.pending,1)
    flist = sprintf('%s\n%s', flist, opts.efunc.pending{i,1});
  end
  error('%s: cannot close equations, function(s) pending:\n%s', ...
    mfilename, flist);
end

%% process input arguments

efunc        = opts.efunc;
cont_pars    = efunc.cont_pars;
tmp          = efunc.par_arrays;
cont_par_idx = coco_par2idx(opts, cont_pars, 'sloppy');

%% augment list of parameters for screen output

efunc.op_idx = [ efunc.op_idx cont_par_idx(tmp.cpbegin:end) ];
efunc.op_idx(efunc.op_idx==0) = [];

% list of all parameters for inclusion in bd
bdp_idx = [efunc.inactive_pars efunc.active_pars efunc.internal_pars ...
  efunc.inequality_pars efunc.regular_pars efunc.singular_pars];
efunc.bdp_idx = sort(bdp_idx);

%% create event structure efunc.ev

ev.pidx       = [];
ev.midx       = [];
ev.par_type   = {};
ev.vals       = [];
ev.point_type = {};
ev.evgroup    = {};
ev.idx        = [];
ev.evsign     = '';
ev.BP_idx     = [];
ev.MX_idx     = [];
ev.SP_idx     = [];

% create copy of indices of remaining parameters
pararrays = { 'regular_pars' 'singular_pars' };
for i = 1:numel(pararrays)
	tmp.(pararrays{i}) = efunc.(pararrays{i});
end

% we ignore events in inactive parameters since these remain constant,
% but include primary continuation parameters
pararrays = { 'primary_pars' 'active_pars' 'inactive_pars' ...
  'internal_pars' 'inequality_pars' 'regular_pars' 'singular_pars' };
partypes  = { 'continuation' 'continuation' 'continuation'  ...
	'continuation' 'continuation' 'regular' 'singular' };

if isfield(efunc, 'events')
  events = efunc.events;
else
	events = [];
end

evnum = 0;

for i=1:numel(events)
	
	% look for parameters name in idx2par
	pars = events(i).par;
	pidx = [];
	for j=1:numel(pars)
		idx = find( strcmp(pars{j}, efunc.idx2par), 1 );
		if isempty(idx)
			error('%s: parameter ''%s'' not found, cannot add events', ...
				mfilename, pars{j});
		end
		pidx = [ pidx ; idx ]; %#ok<AGROW>
	end
	
	% compute type of event parameters
	ptype = {};
	for j=1:numel(pidx)
		type = {};
		for k=1:numel(pararrays)
			if ~isempty(find(tmp.(pararrays{k})==pidx(j),1))
				type = partypes{k};
				break
			end
		end
		if isempty(type)
			ptype = {};
			break;
		else
			ptype = [ ptype type ]; %#ok<AGROW>
		end
	end
	if isempty(ptype)
		continue
	end

	vals       = events(i).vals;
	vnum       = numel(vals);
	evidx      = evnum+(1:vnum);
	point_type = events(i).name;
	evgroup    = {[]};
	evsign     = events(i).sign;
  idx        = i;

	% expand arrays to match size of vals
	o            = ones(vnum, 1);
	point_type   = { point_type{o} };
	if numel(pidx)==1
		pidx       =   pidx(o);
		ptype      = { ptype{o}   };
		evsign     =   evsign(o)';
	else
		evgroup    = { evidx       };
	end
	evgroup      = { evgroup{o}  };
	idx          =   idx(o);

	% update entries in event structure
	ev.pidx       = [ ev.pidx         ; pidx                  ];
  ev.midx       = [ ev.midx         ; efunc.pidx2midx(pidx) ];
	ev.par_type   = [ ev.par_type       ptype                 ];
	ev.vals       = [ ev.vals         ; vals                  ];
	ev.point_type = [ ev.point_type     point_type            ];
	ev.evgroup    = [ ev.evgroup        evgroup               ];
	ev.evsign     = [ ev.evsign         evsign                ];
	ev.idx        = [ ev.idx          ; idx                   ];
	
	evlist      = events(i).evlist;
	ev.(evlist) = [ev.(evlist) evidx];

	evnum       = evnum + vnum;
end

efunc.ev   = ev;
opts.efunc = efunc;

opts = efunc_minit(opts);

end
