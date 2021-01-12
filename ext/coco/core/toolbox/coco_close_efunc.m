function [opts, chart] = coco_close_efunc(opts, cont_pars)
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
% $Id: coco_close_efunc.m 3178 2020-02-20 17:07:40Z hdankowicz $

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

%% process input arguments

if ischar(cont_pars)
  cont_pars = { cont_pars };
end

%% add pending functions

opts = coco_add_pending(opts, 'efunc');
opts = coco_add_pending(opts, 'coco_adjoint');

efunc           = opts.efunc;
efunc.cont_pars = cont_pars;
cont_par_idx    = coco_par2idx(opts, cont_pars, 'sloppy');
miss_par_idx    = find(cont_par_idx==0, 1);

%% compute number of parameters to activate
% ppnum = m + d - n, m+d=efunc.f_dim, n = efunc.x_dim+efunc.p_dim
ppnum = efunc.f_dim-(efunc.x_dim+efunc.p_dim);
if ppnum<0
	% we have more variables than equations
  errmsg = sprintf('%s: cannot close equations, system underdetermined\n', mfilename);
  errmsg = sprintf('%snumber of equations                   : %d\n', errmsg, efunc.f_dim);
  errmsg = sprintf('%snumber of variables+active parameters : %d', errmsg, efunc.x_dim+efunc.p_dim);
	error(errmsg); %#ok<SPERR>
end

%% create copy of indices of potential continuation parameters
pararrays = { 'inactive_pars' 'active_pars' 'internal_pars' 'inequality_pars'};

for i = 1:numel(pararrays)
	tmp.(pararrays{i}) = efunc.(pararrays{i});
end

%% exchange parameters as defined with coco_xchg_pars
if isfield(efunc, 'xchg')
	xchg = coco_par2idx(opts, efunc.xchg);
	tmp  = exchange_pars(tmp, pararrays, xchg);
end

%% exchange internal parameters
pnum  = numel(cont_pars);         % number of parameters passed as argument
ipnum = numel(tmp.internal_pars); % number of internal pars
ipend = min(ipnum, pnum-ppnum);   % exchange internal parameters 1:ipend

if ipend>0 % more cont_pars than required -> exchange internal parameters
  if miss_par_idx<=ppnum+ipend
    error('%s: parameter ''%s'' not found', mfilename, cont_pars{miss_par_idx});
  end
	xchg = [];
	for i = 1:ipend
		idx1 = tmp.internal_pars(i);
		idx2 = cont_par_idx(ppnum+i);
		xchg = [ xchg ; idx1 idx2 ]; %#ok<AGROW>
	end
	tmp = exchange_pars(tmp, pararrays, xchg);
end

%% compute subset of inactive parameters of argument cont_pars
inactive_cont_par_idx = [];
for idx = cont_par_idx
	if any(tmp.inactive_pars == idx) || idx==0 % always treat missing pars inactive
    inactive_cont_par_idx = [ inactive_cont_par_idx idx ]; %#ok<AGROW>
	end
end
cpnum = numel(inactive_cont_par_idx); % number of inactive cont_pars

if cpnum<ppnum
	% we have more equations than variables+parameters
  errmsg = sprintf('%s: cannot close equations, too few parameters activated?\n', mfilename);
  errmsg = sprintf('%snumber of parameters to activate   : %d\n', errmsg, ppnum);
  errmsg = sprintf('%snumber of inactive parameters given: %d', errmsg, cpnum);
	error(errmsg); %#ok<SPERR>
end

%% activate ppnum inactive parameters

tmp.primary_pars = inactive_cont_par_idx(1:ppnum);
if any(tmp.primary_pars==0)
  error('%s: parameter ''%s'' not found', mfilename, cont_pars{miss_par_idx});
end
% for i = 1:ppnum
% 	% verify that cont_par_idx(1:ppnum) are now inactive
% 	if isempty(find(tmp.inactive_pars == cont_par_idx(i), 1))
% 		error('%s: too few inactive parameters specified for continuation', ...
% 			mfilename);
% 	end
% 	add parameter to list primary_pars
% 	tmp.primary_pars = [ tmp.primary_pars inactive_cont_par_idx(i)];
% end

acp_idx          = [ tmp.primary_pars tmp.internal_pars tmp.active_pars ...
  tmp.inequality_pars ];
% bug: acp_idx needs to be reordered according to input list in
% cont_par_idx
% efunc.acp_idx    = acp_idx(arrayfun(@(x) find(x==cont_par_idx,1), acp_idx));
efunc.acp_idx    = acp_idx;
efunc.acp_f_idx  = efunc.pidx2fidx(efunc.acp_idx);

%% check for duplicate active continuation parameters

for i=1:numel(efunc.acp_idx)
	idx = find(efunc.acp_idx == efunc.acp_idx(i));
	if numel(idx) ~= 1
		error('%s: duplicate continuation parameter ''%s''', ...
			mfilename, efunc.idx2par{efunc.acp_idx(i)});
	end
end

%% create initial part of list of parameters for screen output

efunc.op_idx     = [ cont_par_idx(1:ppnum) tmp.internal_pars(ipend+1:end) ];
tmp.cpbegin      = ppnum+1;
efunc.par_arrays = tmp;

%% create permutation vector for xp

efunc.p_idx  = efunc.x_dim+(1:efunc.p_dim+ppnum);
efunc.xp_idx = [efunc.x_idx efunc.p_idx];
efunc.xp_dim = efunc.x_dim+efunc.p_dim+ppnum;

opts.efunc = efunc;

[opts, chart] = efunc_init(opts);

end

%% local functions

function tmp = exchange_pars(tmp, pararrays, xchg)

for i=1:size(xchg,1)

	for j=1:numel(pararrays)
		old_pararray = pararrays{j};
		old_pidx = find(tmp.(old_pararray)==xchg(i,1));
		if ~isempty(old_pidx); break; end
	end
	
	for j=1:numel(pararrays)
		new_pararray = pararrays{j};
		new_pidx = find(tmp.(new_pararray)==xchg(i,2));
		if ~isempty(new_pidx); break; end
	end
	
	pidx = tmp.(old_pararray)(old_pidx);
	tmp.(old_pararray)(old_pidx) = tmp.(new_pararray)(new_pidx);
	tmp.(new_pararray)(new_pidx) = pidx;
  
end

end
