function prob = bvp_sol2seg(prob, oid, varargin)
%BVP_SOL2SEG   Append 'bvp' instance constructed from saved data.
%
% Support restarting continuation from a previously obtained solution,
% stored to disk.
%
% Identical to bvp_v1.
%
% PROB     = BVP_SOL2SEG(PROB, OID, VARARGIN)
% VARARGIN = { RUN [SOID] LAB }
%
% PROB - Continuation problem structure.
% OID  - Target object instance identifier (string).
% RUN  - Run identifier (string).
% SOID - Source object instance identifier (string).
% LAB  - Solution label (integer).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: bvp_sol2seg.m 2839 2015-03-05 17:09:01Z fschild $

ttbid = coco_get_id(oid, 'bvp'); % Create toolbox instance identifier
str  = coco_stream(varargin{:}); % Convert varargin to stream of tokens for argument parsing
run  = str.get;
if ischar(str.peek)
  stbid = coco_get_id(str.get, 'bvp');
else
  stbid = ttbid;
end
lab = str.get;

data = coco_read_solution(stbid, run, lab);
toid = coco_get_id(ttbid, 'seg'); % Target 'coll' object instance identifier
soid = coco_get_id(stbid, 'seg'); % Source 'coll' object instance identifier
prob = coll_sol2seg(prob, toid, run, soid, lab); % Construct 'coll' instance
data = bvp_init_data(prob, ttbid, data); % Build toolbox data
prob = bvp_close_seg(prob, ttbid, data); % Append boundary conditions

end
