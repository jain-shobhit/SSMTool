function prob = ode_isol2hspo(prob, oid, varargin)
%ODE_ISOL2HSPO   Start continuation of multisegment periodic orbits from initial guess.
%
% PROB     = ODE_ISOL2HSPO(PROB, OID, VARARGIN)
% VARARGIN = { FUNCS SIG {T0...} {X0...} [PNAMES] P0 [OPTS] }
% FUNCS    = FUN [ DX [ DP [DXDX [DXDP [DPDP]]]] ]
% FUN      = { @F @E @R }
% DX       = { (@DFDX | '[]') (@DEDX | '[]') (@DRDX | '[]') }
% DP       = { (@DFDP | '[]') (@DEDP | '[]') (@DRDP | '[]') }
% DXDX     = { (@DFDXDX | '[]') (@DEDXDX | '[]') (@DRDXDX | '[]') }
% DXDP     = { (@DFDXDP | '[]') (@DEDXDP | '[]') (@DRDXDP | '[]') }
% DPDP     = { (@DFDPDP | '[]') (@DEDPDP | '[]') (@DRDPDP | '[]') }
% SIG      = { MODE... } { EVENT... } { RESET... }
% OPTS     = { '-hspo-end' | '-end-hspo' | '-var' VECS }
%
% Start a continuation of multisegment periodic orbits, where
%
%   * the i-th segment is governed by an autonomous evolution equation of
%     the form x'=f(x,p;modes{i}) for some non-linear function f,
%   * the i-th segment terminates on the zero-level surface
%     h(x,p;events{i}) for some non-linear function h, and
%   * the initial end point x0 on the i-th segment coincides with the image
%     g(x1,p;reset{i}) of the final end point x1 on the preceding segment
%     for some non-linear function g.
%
% Collectively, the modes, events, and resets arrays are the signature of
% the periodic orbit. Each such array is populated with string labels drawn
% from a finite list.
%
% If the list of possible mode labels consists of at least two elements,
% then the function f should be encoded in a separate m-file using a switch
% statement to distinguish between the vector fields for different mode
% labels. Similar statements apply to the functions h and g for different
% events and reset labels. As an example, the vectorized encoding
%
%   function y = stickslip(x, p, mode)
%
%   switch mode
%     case 'stick'
%       F = p(1,:).^2.*sin(x(3,:)).^2./(1-x(1,:)).^2-c.*x(2,:)-x(1,:);
%
%       y(1,:) = x(2,:);
%       y(2,:) = F;
%       y(3,:) = p(4,:);
%     case 'slip'
%       F = p(1,:).^2.*sin(x(4,:)).^2./(1-x(2,:)).^2-p(2,:).*x(3,:)-x(2,:);
%
%       y(1,:) = -p(3,:)-F/5;
%       y(2,:) = x(3,:);
%       y(3,:) = p(3,:)+6*F/5;
%       y(4,:) = p(4,:);
%   end
%
%   end
%
% defines two distinct vector fields, of different state-space dimension,
% associated with the mode labels 'stick' and 'slip', respectively.
% Similarly, the function
%
%   function y = stickslip_events(x, p, event)
%
%   switch event % corrected typo in Recipes for Continuation, 1st edition, page 240
%     case 'collision'
%       y = 0.5-x(1);
%     case 'phase'
%       y = pi/2-x(3);
%     case 'minsep'
%       y = x(2);
%     case 'rest'
%       y = x(1);
%   end
%
%   end
%
% encodes four distinct event functions associated with the event labels
% 'collision', 'phase', 'minsep', and 'rest', respectively. Finally, the
% function
%
%   function y = stickslip_resets(x, p, reset)
%
%   switch reset
%     case 'bounce'
%       y(1) = (1+p(5))/6*x(2);
%       y(2) = x(1);
%       y(3) = -p(5)*x(2);
%       y(4) = x(3);
%     case 'phase'
%       y(1) = x(1);
%       y(2) = x(2);
%       y(3) = x(3)-pi;
%     case 'turn'
%       y = x;
%     case 'stick'
%       y = x(2:4);
%   end
%
% end
%
% encodes four distinct reset functions associated with the event labels
% 'bounce', 'phase', 'turn', and 'stick', respectively.
%
% On input:
%
% PROB   : Continuation problem structure. This structure is either
%          automatically created and passed on by COCO when using the PO
%          toolbox implicitly, or created with COCO_PROB prior to calling
%          ODE_ISOL2HSPO explicitly. See below for implicit and explicit
%          use of the PO toolbox.
%
%          The continuation problem structure may contain settings defined
%          with COCO_SET, which will influence the behavior of the PO
%          toolbox and each 'coll' instance embedded in an 'hspo' instance.
%          Execute HSPO_SETTINGS on the Matlab command line to see an
%          overview of PO toolbox settings specific to multisegment
%          periodic orbits.
%
% OID    : Object instance identifier (string). The corresponding toolbox
%          instance identifier is coco_get_id(OID, 'hspo'). Pass the empty
%          string '' for a simple continuation of multisegment periodic
%          orbits. Pass a non-trivial object identifier if an instance of
%          the PO toolbox is part of a composite continuation problem.
%
% F      : Function handle to vector field.
% E      : Function handle to event function.
% R      : Function handle to reset function.
% DFDX   : Function handle to Jacobian of vector field w.r.t. state
%          variables (optional, may be empty).
% DEDX   : Function handle to Jacobian of event function w.r.t. state
%          variables (optional, may be empty).
% DRDX   : Function handle to Jacobian of reset function w.r.t. state
%          variables (optional, may be empty).
% DFDP   : Function handle to Jacobian of vector field w.r.t. problem
%          parameters (optional, may be empty).
% DEDP   : Function handle to Jacobian of event function w.r.t. problem
%          parameters (optional, may be empty).
% DRDP   : Function handle to Jacobian of reset function w.r.t. problem
%          parameters (optional, may be empty).
% MODE   : Mode identifier (string).
% EVENT  : Event identifier (string).
% RESET  : Reset identifier (string).
% T0     : Segment-by-segment initial time discretization (cell array).
% X0     : Segment-by-segment initial sampled time history for state
%          variables (cell array).
% PNAMES : Cell array of string labels for continuation parameters
%          (optional). These string labels will be mapped onto the problem
%          parameters as lab{1}->p(1), lab{2}->p(2), ... If string labels
%          are given, their number must match the number of problem
%          parameters. If only one problem parameter is present, a string
%          may be given instead of a cell array of strings. Omitting this
%          argument is only relevant for embedded problems. In general,
%          this argument must be given.
% P0     : Initial solution guess for problem parameters. Pass an empty
%          array to indicate absence of problem parameters.
%
% OPTS   : '-hspo-end', '-end-hspo', and '-var' VECS (optional, multiple
%          options may be given). Either '-hspo-end' or '-end-hspo' marks
%          the end of input to ODE_ISOL2HSPO. The option '-var' indicates
%          the inclusion of the variational problem for each trajectory
%          segment, where the initial solution guess for the perturbations
%          to the initial conditions  of the i-th trajectory segment is
%          given by the content of the i-th element of the cell array VECS.
%
% Implicit use of PO toolbox
% --------------------------
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'hspo', HSPO_ARGS, CONT_ARGS);
%
% will call the function ODE_ISOL2HSPO implicitly, where the arguments
% HSPO_ARGS following the string 'coll' will be passed to ODE_ISOL2HSPO.
% Any remaining arguments CONT_ARGS will be passed to the continuation
% method. When using the default 1-dimensional continuation method,
% CONT_ARGS defines the list of continuation parameters and computational
% domains.
%
% A call to COCO of the form
%
%   coco(prob, RUN, 'ode', 'isol', 'hspo', HSPO_ARGS, CONT_ARGS);
%
% differs from the above in that the existing continuation problem
% structure prob is used in the call to ODE_ISOL2HSPO. Non-default values
% for HSPO and COLL toolbox settings may be assigned to prob before this
% call.
%
% Explicit use of PO toolbox
% --------------------------
% A calling sequence like
%
%   prob = coco_prob();
%   ...
%   prob = ode_isol2hspo(prob, '', HSPO_ARGS);
%   ...
%   coco(prob, RUN, [], CONT_ARGS);
%
% uses the function ODE_ISOL2HSPO explicitly to add a 'hspo' instance to
% the continuation problem structure prob. This calling sequence is
% equivalent to the implicit use of the PO toolbox, but allows for
% inclusion of user-defined monitor functions between the calls to
% ODE_ISOL2HSPO and COCO, as demonstrated in the demos.
%
% See also: COCO_SET, HSPO_SETTINGS, HSPO_READ_SOLUTION, HSPO_ADD

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_isol2segs.m 2839 2015-03-05 17:09:01Z fschild $

grammar   = 'F [DFDX [DFDP [DFDXDX [DFDXDP [DFDPDP]]]]] MODES EVENTS RESETS T0 X0 [PNAMES] P0 [OPTS]';
args_spec = {
  'F', 'cell',     '{@}',      'fhan',         [], 'read', {}
  'DFDX', 'cell',  '{@|[]}',   'dfdxhan', {[] [] []}, 'read', {}
  'DFDP', 'cell',  '{@|[]}',   'dfdphan', {[] [] []}, 'read', {}
  'DFDXDX', 'cell',  '{@|[]}', 'dfdxdxhan', {[] [] []}, 'read', {}
  'DFDXDP', 'cell',  '{@|[]}', 'dfdxdphan', {[] [] []}, 'read', {}
  'DFDPDP', 'cell',  '{@|[]}', 'dfdpdphan', {[] [] []}, 'read', {}
  'MODES', 'cell',   '{str}',     'modes',         [], 'read', {}
  'EVENTS', 'cell',   '{str}',    'events',         [], 'read', {}
  'RESETS', 'cell',   '{str}',    'resets',         [], 'read', {}
  'T0', 'cell', '{[num]}',        't0',         [], 'read', {}
  'X0', 'cell', '{[num]}',        'x0',         [], 'read', {}
  'PNAMES', 'cell',   '{str}',    'pnames',         {}, 'read', {}
  'P0',     '',   '[num]',        'p0',         [], 'read', {}
  };
opts_spec = {
  '-hspo-end',     '', '',  'end', {}
  '-end-hspo',     '', '',  'end', {}
  '-var', 'vecs', {}, 'read', {}
  };
[args, opts] = coco_parse(grammar, args_spec, opts_spec, varargin{:});

tbid = coco_get_id(oid, 'hspo');

nos = [numel(args.fhan) numel(args.dfdxhan) numel(args.dfdphan) ...
  numel(args.dfdxdxhan) numel(args.dfdxdphan) numel(args.dfdpdphan)];
assert(all(nos==3), '%s: incomplete function handle specification', tbid);
nos = [numel(args.modes) numel(args.events) numel(args.resets) ...
  numel(args.t0) numel(args.x0)];
assert(~any(diff(nos)), '%s incompatible segment specification', tbid);
assert(isempty(args.pnames) || numel(args.p0)==numel(args.pnames), ...
  '%s: incompatible number of elements for ''p0'' and ''pnames''', tbid);

nsegs = nos(1);
coll  = cell(1,nsegs);
if ~isempty(opts.vecs)
  assert(numel(opts.vecs)==nsegs, ...
    '%s: incompatible specification of arrays of perturbations', ...
    mfilename);
else
  opts.vecs = cell(1,nsegs);
end
for i=1:nsegs
  mode      = args.modes{i};
  fhan      = cfn(args.fhan{1}, mode);
  dfdxhan   = cfn(args.dfdxhan{1}, mode);
  dfdphan   = cfn(args.dfdphan{1}, mode);
  dfdxdxhan = cfn(args.dfdxdxhan{1}, mode);
  dfdxdphan = cfn(args.dfdxdphan{1}, mode);
  dfdpdphan = cfn(args.dfdpdphan{1}, mode);
  coll{i}   = {fhan, dfdxhan, dfdphan, dfdxdxhan, dfdxdphan, dfdpdphan, ...
    args.t0{i}, args.x0{i}, args.p0, '-var', opts.vecs{i}};
end

data = hspo_init_data(prob, args, oid); % initialize 'hspo' instance data structure
data.hspo_orb = init_data(args, oid);
if data.hspo.bifus
  prob = coco_set(prob, coco_get_id(oid, 'coll'), 'var', true);
end
prob = ode_isol2bvp(prob, coco_get_id(tbid, 'orb'), coll, ...
  args.pnames, @bc, @bc_DFDX, @bc_DFDXDX, data);
prob = hspo_add(prob, data);
prob = ode_add_tb_info(prob, oid, 'hspo', 'hspo', 'hspo', hspo_sol_info());

end

function fhan = cfn(fhan, mode)
%cfn   Create 'coll'-compatible function handle.
%
% Wrapper to accommodate that 'coll' compatible vector fields depend only
% on x and p. Different vector fields are identified by the mode label
% and handled within the functional encoding using a switch statement.

if isa(fhan, 'function_handle'), fhan = @(x,p) fhan(x,p,mode); end
end

function data = init_data(args, oid)
%INIT_DATA   Initialize boundary conditions data for an instance of 'hspo'.

% On output:
%
% The DATA output argument contains the fields
%
%    nsegs  : number of trajectory segments
%    x0_idx : cell array with context-independent index arrays of initial
%             end points of trajectory segments in vector x0 of
%             continuation variables in boundary conditions
%    x1_idx : cell array with context-independent index arrays of final
%             end points of trajectory segments in vector x1 of
%             continuation variables in boundary conditions
%    dim    : array of state-space dimensions
%    cdim   : total state-space dimension
%    pdim   : number of problem parameters
%    Jnz, Jcnz, Jrows, and Jcols : index arrays used to allocate space in
%             Jacobian of boundary conditions

nsegs       = numel(args.modes);
data.nsegs  = nsegs;
data.x0_idx = cell(1,nsegs);
data.x1_idx = cell(1,nsegs);
cdim        = 0;
dim         = zeros(1,nsegs);
for i=1:nsegs
  dim(i)         = size(args.x0{i},2);
  data.x0_idx{i} = cdim + (1:dim(i))';
  data.x1_idx{i} = cdim + (1:dim(i))';
  cdim           = cdim + dim(i);
end
data.dim  = dim;
data.adim = [0 cumsum(dim)];
data.cdim = cdim;  % Total state-space dimension.
pdim      = numel(args.p0);
data.pdim = pdim;  % Number of problem parameters
data.m    = nsegs+cdim;
data.n    = nsegs+2*cdim+pdim;
Jrows = [];
Jcols = [];
Jnz   = zeros(1,nsegs);
off  = 0;
for i=1:nsegs
  j      = mod(i,nsegs)+1;
  Jrows  = [Jrows; repmat(off+1, [dim(i)+pdim 1])]; %#ok<AGROW>
  Jrows  = [Jrows; repmat(off+1+(1:dim(j))', [1+dim(i)+pdim 1])]; %#ok<AGROW>
  Jcols  = [Jcols; nsegs+cdim+data.x1_idx{i}; nsegs+2*cdim+(1:pdim)']; %#ok<AGROW>
  c2     = repmat(nsegs+cdim+data.x1_idx{i}', [dim(j) 1]);
  c3     = repmat(nsegs+2*cdim+(1:pdim), [dim(j) 1]);
  Jcols  = [Jcols; nsegs+data.x0_idx{j}; c2(:); c3(:)]; %#ok<AGROW>
  off    = off + dim(j) + 1;
  Jnz(i) = (dim(j)+1)*(dim(i)+pdim)+dim(j);
end
data.Jnz   = Jnz;
data.Jcnz  = sum(Jnz);
data.Jrows = Jrows; % Row index array for sparse Jacobian
data.Jcols = Jcols; % Column index array for sparse Jacobian

dJrows = [];
dJcols = [];
dJnz   = zeros(1,nsegs);
off   = 0;
for i=1:nsegs
  j = mod(i,nsegs)+1;
  dJrows = [dJrows; repmat(off+1, [(dim(i)+pdim)^2 1])]; %#ok<AGROW>
  dJrows = [dJrows; repmat(off+1+(1:dim(j))', [(dim(i)+pdim)^2 1])]; %#ok<AGROW>
  % dhdxdx
  cols   = repmat(nsegs+cdim+data.x1_idx{i}, [1 dim(i)]);
  cols   = cols + ...
    repmat(data.n*(nsegs+cdim+data.adim(i)):data.n: ...
    data.n*(nsegs+cdim+data.adim(i+1)-1), [dim(i) 1]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  % dhdxdp
  cols   = repmat(nsegs+cdim+data.x1_idx{i}, [1 pdim]);
  cols   = cols + ...
    repmat(data.n*(nsegs+2*cdim):data.n: ...
    data.n*(nsegs+2*cdim+pdim-1), [dim(i) 1]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  % dhdpdx
  cols   = repmat(nsegs+cdim+data.x1_idx{i}, [1 pdim]);
  cols   = cols + ...
    repmat(data.n*(nsegs+2*cdim):data.n: ...
    data.n*(nsegs+2*cdim+pdim-1), [dim(i) 1]);
  cols   = cols';
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  % dhdpdp
  cols   = repmat(nsegs+2*cdim+(1:pdim)', [1 pdim]);
  cols   = cols + ...
    repmat(data.n*(nsegs+2*cdim):data.n: ...
    data.n*(nsegs+2*cdim+pdim-1), [pdim 1]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  
  % dgdxdx
  cols   = repmat(nsegs+cdim+data.x1_idx{i}', [dim(j) 1 dim(i)]);
  cols   = cols + ...
    repmat(reshape(data.n*(nsegs+cdim+data.adim(i)):data.n: ...
    data.n*(nsegs+cdim+data.adim(i+1)-1), [1 1 dim(i)]), [dim(j) dim(i) 1]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  % dgdxdp
  cols   = repmat(nsegs+cdim+data.x1_idx{i}', [dim(j) 1 pdim]);
  cols   = cols + ...
    repmat(reshape(data.n*(nsegs+2*cdim):data.n: ...
    data.n*(nsegs+2*cdim+pdim-1), [1 1 pdim]), [dim(j) dim(i) 1]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  % dgdpdx
  cols   = repmat(nsegs+cdim+data.x1_idx{i}', [dim(j) 1 pdim]);
  cols   = cols + ...
    repmat(reshape(data.n*(nsegs+2*cdim):data.n: ...
    data.n*(nsegs+2*cdim+pdim-1), [1 1 pdim]), [dim(j) dim(i) 1]);
  cols   = permute(cols, [1 3 2]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  % dgdpdp
  cols   = repmat(nsegs+2*cdim+(1:pdim), [dim(j) 1 pdim]);
  cols   = cols + ...
    repmat(reshape(data.n*(nsegs+2*cdim):data.n: ...
    data.n*(nsegs+2*cdim+pdim-1), [1 1 pdim]), [dim(j) pdim 1]);
  dJcols = [dJcols; cols(:)]; %#ok<AGROW>
  
  off = off+dim(j)+1;
  dJnz(i) = (1+dim(j))*(dim(i)+pdim)^2;
end
data.dJnz   = dJnz;
data.dJcnz  = sum(dJnz);
data.dJrows = dJrows; % Row index array for sparse Jacobian
data.dJcols = dJcols; % Column index array for sparse Jacobian

data.fid  = coco_get_id(oid, 'hspo');

end

function y = bc(data, T, x0, x1, p) %#ok<INUSL>
%BC   Multi-segment periodic boundary conditions.

pr = data.pr;
hspo_orb = pr.hspo_orb;
nsegs    = hspo_orb.nsegs;

y = zeros(hspo_orb.cdim+nsegs,1);
off = 0;
for i=1:nsegs
  pt0 = x0(hspo_orb.x0_idx{mod(i,nsegs)+1});
  pt1 = x1(hspo_orb.x1_idx{i});
  y(off + (1:hspo_orb.dim(mod(i,nsegs)+1) + 1)) = ...
    [pr.event_F(pr, pt1, p, pr.events{i}); ...
    pt0 - pr.reset_F(pr, pt1, p, pr.resets{i})];
  off = off + hspo_orb.dim(mod(i,nsegs)+1) + 1;
end

end

function J = bc_DFDX(data, T, x0, x1, p) %#ok<INUSL>
%BC_DFDX   Linearization of multi-segment periodic boundary conditions.

pr = data.pr;
hspo_orb = pr.hspo_orb;
nsegs = hspo_orb.nsegs;

vals = zeros(1,hspo_orb.Jcnz);
off  = 0;
for i=1:nsegs
  pt1  = x1(hspo_orb.x1_idx{i});
  dhdx = pr.event_DFDX(pr, pt1, p, pr.events{i});
  dhdp = pr.event_DFDP(pr, pt1, p, pr.events{i});
  dgdx = pr.reset_DFDX(pr, pt1, p, pr.resets{i});
  dgdp = pr.reset_DFDP(pr, pt1, p, pr.resets{i});
  dim  = hspo_orb.dim(mod(i,nsegs)+1);
  vals(off + (1:hspo_orb.Jnz(i))) = ...
    [dhdx(:); dhdp(:); ones(dim,1); -dgdx(:); -dgdp(:)];
  off  = off + hspo_orb.Jnz(i);
end
J = sparse(hspo_orb.Jrows, hspo_orb.Jcols, vals, hspo_orb.m, hspo_orb.n);

end

function dJ = bc_DFDXDX(data, T, x0, x1, p) %#ok<INUSL>
%BC_DFDXDX   2nd linearization of multi-segment periodic boundary conditions.

pr = data.pr;
hspo_orb = pr.hspo_orb;
nsegs = hspo_orb.nsegs;

vals = zeros(1,hspo_orb.dJcnz);
off  = 0;
for i=1:nsegs
  pt1  = x1(hspo_orb.x1_idx{i});
  dhdxdx = pr.event_DFDXDX(pr, pt1, p, pr.events{i});
  dhdxdp = pr.event_DFDXDP(pr, pt1, p, pr.events{i});
  dhdpdp = pr.event_DFDPDP(pr, pt1, p, pr.events{i});
  dgdxdx = pr.reset_DFDXDX(pr, pt1, p, pr.resets{i});
  dgdxdp = pr.reset_DFDXDP(pr, pt1, p, pr.resets{i});
  dgdpdp = pr.reset_DFDPDP(pr, pt1, p, pr.resets{i});
  vals(off + (1:hspo_orb.dJnz(i))) = ...
    [dhdxdx(:); dhdxdp(:); dhdxdp(:); dhdpdp(:); -dgdxdx(:); ...
    -dgdxdp(:); -dgdxdp(:); -dgdpdp(:)];
  off  = off + hspo_orb.dJnz(i);
end
dJ = sparse(hspo_orb.dJrows, hspo_orb.dJcols, vals, hspo_orb.m, hspo_orb.n^2);

end
