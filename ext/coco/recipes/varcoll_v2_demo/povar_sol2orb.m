function [prob data sol] = povar_sol2orb(prob, oid, varargin)
%POVAR_SOL2ORB   Append 'po' instance and 'varcoll' instance constructed from stored data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% [PROB DATA SOL] = POVAR_SOL2ORB(PROB, OID, VARARGIN)
% VARARGIN        = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% DATA - 'varcoll' instance data.
% SOL  - Fundamental solution to variational equations.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: povar_sol2orb.m 2839 2015-03-05 17:09:01Z fschild $

str   = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run   = str.get;
if ischar(str.peek)
  soid = str.get;
else
  soid = oid;
end
lab = str.get;

stbid = coco_get_id(soid, 'po.seg.var');
[data sol] = coco_read_solution(stbid, run, lab); % Extract 'varcoll' data and fundamental solution
prob = po_sol2orb(prob, oid, run, soid, lab); % Reconstruct 'po' instance
toid = coco_get_id(oid, 'po.seg');
prob = var_coll_add(prob, toid, data.dfdxdxhan, data.dfdxdphan); % Append variational zero problem

end
