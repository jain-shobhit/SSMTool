function [sols, data] = tor_read_solution(oid, run, lab)
%TOR_READ_SOLUTION   Read 'tor' solution and toolbox data from disk.
%
% Extract data and chart structures associated with 'tor' toolbox instance
% identifier from solution file and construct solution structure.
%
% [SOL DATA] = TOR_READ_SOLUTION(OID, RUN, LAB)
%
% SOL  - Solution (struct).
% DATA - Toolbox data (struct).
% OID  - Object instance identifier (string).
% RUN  - Run identifier (string).
% LAB  - Solution label (integer).

tbid   = coco_get_id(oid, 'tor');
[sol, data] = bvp_read_solution(tbid, run, lab);
N  = data.nsegs;
% set unified time mesh for all segments
tbp = sol{1}.tbp;
assert(abs(tbp(1))<1e-4,'initial time is not zero');
flags = false(N,1);
for j=1:N
    n  = mod(j-1,N)+1;
    tj = sol{n}.tbp;
    assert(abs(tj(1))<1e-4,'initial time is not zero');
    assert(abs(tj(end)-tbp(end))<1e-3*abs(tj(end)),'final time is not matched in segments of torus');
    if numel(tj)==numel(tbp)
        flags(j) = norm(tj-tbp)<1e-3*norm(tbp);
    end
    if numel(tj)>numel(tbp)
        tbp = tj; % pick the time mesh with most points
    end
end

sameMesh = all(flags);

% interpolation of trajectories at the unified mesh
nt  = numel(tbp);
xbp = zeros(nt, data.bc_data.dim, N+1);
for j=1:N+1
    n    = mod(j-1,N)+1; % the last one is the same as the first one for torus plotting
    tbpj = sol{n}.tbp;
    xbpj = sol{n}.xbp;
    if ~sameMesh
        xbpj = interp1(tbpj, xbpj, tbp, 'pchip');
    end
    xbp(:,:,j) = xbpj;
end

% write output
sols = struct();
sols.tbp = tbp;
sols.xbp = xbp;
sols.p   = sol{1}.p;

end
