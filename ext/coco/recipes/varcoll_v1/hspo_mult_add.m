function prob = hspo_mult_add(prob, segsoid)
%HSPO_MULT_ADD   Add slot and/or monitor function computing Floquet multipliers of hybrid periodic orbit.
%
% Append slot and/or monitor function to each segment of an 'hspo' instance.
%
% HSPO_MULT_ADD(PROB, SEGSOID)
%
% PROB    - Continuation problem structure.
% SEGSOID - Multi-segment boundary-value problem object instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: hspo_mult_add.m 2839 2015-03-05 17:09:01Z fschild $

cids  = {};
vids  = {};
msid  = coco_get_id(segsoid, 'msbvp'); % 'msbvp' toolbox instance identifier
fdata = coco_get_func_data(prob, msid, 'data'); % Extract 'msbvp' instance data
data  = fdata.bc_data; % 'hspo' instance data
for i=1:data.nsegs
  soid = coco_get_id(msid, sprintf('seg%d', i));
  cids = [cids {coco_get_id(soid, 'coll')}];
  vids = [vids {coco_get_id(soid, 'var')}]; 
end
data.msid   = msid;
data.cids   = cids; % 'coll' object instance identifiers
data.vids   = vids; % 'varcoll' object instance identifiers
data.p_idx  = numel(fdata.x1_idx)+(1:numel(fdata.p_idx));

% This is the solution from Recipes for Continuation, 1st edition, page
% 263, albeit with the corrected implementation of hspo_mult_eigs_bddat,
% which eliminates the incorrect use of function data.

data.mnames = coco_get_id(msid, 'multipliers'); % Column header
prob = coco_add_slot(prob, data.mnames, @hspo_mult_eigs_bddat, data, 'bddat');

% This is an alternative solution that makes correct use of function data.

fdata = coco_get_func_data(prob, data.cids{1}, 'data');
mfid  = coco_get_def_par_names(data.mnames, 1:fdata.dim);
uidx  = coco_get_func_data(prob, msid, 'uidx');
prob  = coco_add_func(prob, data.mnames, @hspo_mult_eigs, data, 'regular', ...
  mfid, 'requires', data.vids, 'uidx', uidx); % Explicit dependence on 'varcoll'

end

function [data y] = hspo_mult_eigs(prob, data, u)
%HSPO_MULT_EIGS   Compute Floquet multipliers.
%
% This function is not included in Recipes for Continuation, 1st edition.
%
% Compute eigenvalues of product of transfer matrices corresponding to the
% Jacobians of the maps from each segment end point at t=0 to next segment
% end point at t=0.

fdata = coco_get_func_data(prob, data.msid, 'data'); % 'msbvp' instance data
P     = hspo_P(prob, data, u([fdata.x1_idx; fdata.p_idx])); % Jacobian of transfer matrix
M = P{1};
for i=2:data.nsegs
  M = P{i}*M;
end

y = eig(M); % Includes trivial eigenvalue at 0

end

function [data res] = hspo_mult_eigs_bddat(prob, data, command, varargin)
%HSPO_MULT_EIGS_BDDAT   Slot function: add Floquet multipliers to bifurcation data.
%
% Compute eigenvalues of product of transfer matrices corresponding to the
% Jacobians of the maps from each segment end point at t=0 to next segment
% end point at t=0.
 
switch command
  case 'init'
    res = {data.mnames};
  case 'data'
    % The following eight lines from Recipes for Continuation, 1st edition,
    % page 262, make incorrect use of function data to extract the
    % Jacobians of the time-T flows. Notably, the available function data
    % of each of the instances of the 'varcoll' toolbox that are accessed
    % in hspo_P need not be associated with the chart that is being flushed
    % to the bifurcation data.

    % [fdata uidx] = coco_get_func_data(prob, data.msid, 'data', 'uidx'); % 'msbvp' instance data and index array
    % chart = varargin{1};
    % u     = chart.x(uidx);
    % P     = hspo_P(prob, data, u([fdata.x1_idx; fdata.p_idx])); % Jacobian of transfer matrix
    % M = P{1};
    % for i=2:data.nsegs
    %   M = P{i}*M;
    % end
    
    % The following 21 lines provide a correct solution to the problem of
    % computing the Jacobians of the Poincare maps, without calling hspo_P.
    % The implementation guarantees a just-in-time solution of each of the
    % variational equations.

    chart = varargin{1};
    M = eye(data.dim(1));
    for i=1:data.nsegs
      [fdata uidx] = coco_get_func_data(prob, data.cids{i}, 'data', 'uidx');
      u  = chart.x(uidx);
      M1 = var_eval_sol(fdata, u);
      x  = u(fdata.x1_idx); % Segment end point at t=1
      p  = u(fdata.p_idx);
      fs = fdata.fhan(x, p);
      if ~isempty(data.dhdxhan)
        hx = data.dhdxhan(x, p, data.events{i});
      else
        hx = coco_ezDFDX('f(x,p)v', data.hhan, x, p, data.events{i});
      end
      if ~isempty(data.dgdxhan)
        gx = data.dgdxhan(x, p, data.resets{i});
      else
        gx = coco_ezDFDX('f(x,p)v', data.ghan, x, p, data.resets{i});
      end
      M  = gx*(eye(fdata.dim)-(fs*hx)/(hx*fs))*M1*M;
    end
    
    res = {eig(M)}; % Includes trivial eigenvalue at 0
end

end

function P = hspo_P(prob, data, u)
%HSPO_P   Compute collection of transfer matrices.
%
% For each segment, extract Jacobian of time-T map and premultiply by
% saltation matrix (correcting for difference in time-of-flight to event
% surface, and the imposition of the reset) to obtain transfer matrix.
%
% This function is obsolete.
%
% P = HSPO_P(PROB, DATA, U)
%
% P    - Cell array of transfer matrices.
% PROB - Continuation problem structure.
% DATA - hspo_mult_eigs_bddat function data.
% U    - Array of end points of segments at t=1 and problem parameters.

P = cell(1,data.nsegs);
p = u(data.p_idx); % Problem parameters
for i=1:data.nsegs
  fdata = coco_get_func_data(prob, data.vids{i}, 'data');
  M1    = fdata.M(fdata.M1_idx,:); % Jacobian of time-T map
  dim   = fdata.dim;
  
  fdata = coco_get_func_data(prob, data.cids{i}, 'data');
  x     = u(data.x1_idx{i}); % Segment end point at t=1
  fs    = fdata.fhan(x, p);
  if ~isempty(data.dhdxhan)
    hx = data.dhdxhan(x, p, data.events{i});
  else
    hx = coco_ezDFDX('f(x,p)v', data.hhan, x, p, data.events{i});
  end
  if ~isempty(data.dgdxhan)
    gx = data.dgdxhan(x, p, data.resets{i});
  else
    gx = coco_ezDFDX('f(x,p)v', data.ghan, x, p, data.resets{i});
  end
  P{i}  = gx*(eye(dim)-(fs*hx)/(hx*fs))*M1; % Multiplication by saltation matrix
end

end
