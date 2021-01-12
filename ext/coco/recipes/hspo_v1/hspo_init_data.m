function data = hspo_init_data(fhan, dfdxhan, dfdphan, modes, events, resets, x0, p0)
%HSPO_INIT_DATA   Initialize boundary conditions data for an instance of 'msbvp'.
%
% Populate remaining fields of the toolbox data structure used by 'hspo'
% function objects.
%
% DATA = HSPO_INIT_DATA(FHAN, DFDXHAN, DFDPHAN, MODES, EVENTS, RESETS, X0, P0)
%
% FHAN    - Cell array of function handles to vector field, event function,
%           and reset function.
% DFDXHAN - Optional cell array of function handles to Jacobians w.r.t.
%           problem variables.
% DFDPHAN - Optional cell array of function handles to Jacobians w.r.t.
%           problem parameters.
% MODES   - Sequence of mode identifiers (cell).
% EVENTS  - Sequence of event identifiers (cell).
% RESETS  - Sequence of reset identifiers (cell).
% X0      - Collection of array of state vectors at mesh points.
% P0      - Initial solution guess for parameter values.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_init_data.m 2839 2015-03-05 17:09:01Z fschild $

data.hhan    = fhan{2};
data.dhdxhan = dfdxhan{2};
data.dhdphan = dfdphan{2};
data.ghan    = fhan{3};
data.dgdxhan = dfdxhan{3};
data.dgdphan = dfdphan{3};
data.modes   = modes;
data.events  = events;
data.resets  = resets;

nsegs       = numel(events);
data.nsegs  = nsegs;         % Number of segments
data.x0_idx = cell(1,nsegs); % Index arrays for trajectory end points at t=0
data.x1_idx = cell(1,nsegs); % Index arrays for trajectory end points at t=1
cdim        = 0;           
dim         = zeros(1,nsegs);
data.dim    = dim; % State-space dimensions
for i=1:nsegs
  dim(i)         = size(x0{i},2);
  data.dim(i)    = dim(i);
  data.x0_idx{i} = cdim+(1:dim(i))';
  data.x1_idx{i} = cdim+(1:dim(i))';
  cdim           = cdim+dim(i);
end
data.cdim = cdim;  % Total state-space dimension.
rows = [];
cols = [];
off  = 0;
pdim = numel(p0);
data.pdim = pdim;  % Number of problem parameters
for i=1:nsegs
  rows = [rows; repmat(off+1, [dim(i)+pdim 1])];
  rows = [rows; repmat(off+1+(1:dim(mod(i,data.nsegs)+1))', ...
    [1+dim(i)+pdim 1])];
  cols = [cols; nsegs+cdim+data.x1_idx{i}; nsegs+2*cdim+(1:pdim)'];
  c2   = repmat(nsegs+cdim+data.x1_idx{i}', ...
    [dim(mod(i,data.nsegs)+1) 1]);
  c3   = repmat(nsegs+2*cdim+(1:pdim), [dim(mod(i,data.nsegs)+1) 1]);
  cols = [cols; nsegs+data.x0_idx{mod(i,data.nsegs)+1}; c2(:); c3(:)];
  off  = off+dim(mod(i,data.nsegs)+1)+1;
end
data.rows = rows; % Row index array for sparse Jacobian
data.cols = cols; % Column index array for sparse Jacobian

end
