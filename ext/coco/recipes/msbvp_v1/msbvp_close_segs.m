function prob = msbvp_close_segs(prob, tbid, data)
%MSBVP_CLOSE_SEGS   Append an instance of 'msbvp' to problem.
%
% Add boundary conditions to multiple instances of 'coll'.
%
% PROB = MSBVP_CLOSE_SEGS(PROB, TBID, DATA)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_close_segs.m 2839 2015-03-05 17:09:01Z fschild $

if ~isempty(data.bc_update) % Optional inclusion of boundary conditions function data
  data.tbid = tbid;
  data = coco_func_data(data); % Convert to func_data class for shared access
  prob = coco_add_slot(prob, tbid, @msbvp_bc_update, data, 'update');
end
T_idx  = zeros(data.nsegs,1);
x0_idx = [];
x1_idx = [];
s_idx  = cell(1, data.nsegs);
for i=1:data.nsegs
  fid      = coco_get_id(tbid,sprintf('seg%d.coll', i)); % Create 'coll' toolbox instance identifier
  [fdata uidx] = coco_get_func_data(prob, fid, 'data', 'uidx'); % Extract 'coll' data structure and context-dependent index array
  T_idx(i) = uidx(fdata.T_idx);            % Subset to T
  x0_idx   = [x0_idx; uidx(fdata.x0_idx)]; % Subset to x0
  x1_idx   = [x1_idx; uidx(fdata.x1_idx)]; % Subset to x1
  s_idx{i} = uidx(fdata.p_idx);            % Subset to p
end
uidx = [T_idx; x0_idx; x1_idx; s_idx{1}]; % Use only one copy of problem parameters
if isempty(data.dfbcdxhan) % Optional inclusion of explicit Jacobian of boundary conditions
  prob = coco_add_func(prob, tbid, @msbvp_F, data, ...
    'zero', 'uidx', uidx);
else
  prob = coco_add_func(prob, tbid, @msbvp_F, @msbvp_DFDU, data, ...
    'zero', 'uidx', uidx);
end
for i=2:data.nsegs % Glue redundant copies of problem parameters
  fid  = coco_get_id(tbid, sprintf('shared%d', i-1));
  prob = coco_add_glue(prob, fid, s_idx{1}, s_idx{i});
end
if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  prob = coco_add_pars(prob, fid, s_idx{1}, data.pnames);
end
prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

function [data y] = msbvp_F(prob, data, u)
%MSBVP_F   COCO-compatible wrapper to boundary conditions.
%
% Expects encoding of boundary conditions as function of T, x0, x1, and p.

T  = u(data.T_idx);  % Extract interval lengths
x0 = u(data.x0_idx); % Extract trajectory end points at t=0
x1 = u(data.x1_idx); % Extract trajectory end points at t=1
p  = u(data.p_idx);  % Extract single copy of problem parameters

y  = data.fbchan(data.bc_data, T, x0, x1, p);

end

function [data J] = msbvp_DFDU(prob, data, u)
%MSBVP_DFDU   COCO-compatible wrapper to linearization of boundary conditions.
%
% Expects encoding of boundary conditions as function of T, x0, x1, and p.

T  = u(data.T_idx);  % Extract interval lengths
x0 = u(data.x0_idx); % Extract trajectory end points at t=0
x1 = u(data.x1_idx); % Extract trajectory end points at t=1
p  = u(data.p_idx);  % Extract single copy of problem parameters

J  = data.dfbcdxhan(data.bc_data, T, x0, x1, p);

end

function data = msbvp_bc_update(prob, data, cseg, varargin)
%MSBVP_BC_UPDATE   COCO-compatible wrapper to boundary condition function data update function.
%
% Use information about current solution to update function data
% parameterizing the execution of the boundary condition zero function.

uidx = coco_get_func_data(prob, data.tbid, 'uidx'); % Context-dependent index set
u    = cseg.src_chart.x(uidx); % Current chart
T    = u(data.T_idx);  % Extract interval lengths
x0   = u(data.x0_idx); % Extract trajectory end points at t=0
x1   = u(data.x1_idx); % Extract trajectory end points at t=1
p    = u(data.p_idx);  % Extract problem parameters
data.bc_data = data.bc_update(data.bc_data, T, x0, x1, p);

end
