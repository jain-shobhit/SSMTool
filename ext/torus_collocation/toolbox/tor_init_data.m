function data = tor_init_data(args,data)
%TOR_INIT_DATA   Initialize toolbox data for an instance of 'tor'.
%
% Populate remaining fields of the toolbox data structure used by 'tor'
% function objects.
%
% DATA = TOR_INIT_DATA(ARGS,DATA)
%
% ARGS - Arguments.
% DATA - Toolbox data structure.

[~,dim,nsegs] = size(args.x0);
assert(mod(nsegs,2)==1, 'the number of segments is not an odd number');
N = (nsegs-1)/2;

% Construct boundary conditions data, Fourier transform and rotation matrix
Th = 2*pi*(0:2*N)/(2*N+1);
Th = kron(1:N, Th');
F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);

data.dim     = dim;
data.N       = N;
data.nsegs   = nsegs;
data.Fs      = F;
data.F       = kron(F, eye(dim));
if data.autonomous
    data.fhan = @(t,x,p) args.fhan(x,p(1:end-3,:)); % with om1, om2 and varrho added as the last three parameters
else
    data.fhan = @(t,x,p) args.fhan(t,x,p(1:end-3,:));
end
if isempty(args.dfdxhan)
    data.dfdxhan = [];
else
    if data.autonomous
        data.dfdxhan = @(t,x,p) args.dfdxhan(x,p(1:end-3,:));
    else
        data.dfdxhan = @(t,x,p) args.dfdxhan(t,x,p(1:end-3,:));
    end
    
end
if isempty(args.dfdphan)
    data.dfdphan = [];
else
    if data.autonomous
        data.dfdphan = @(t,x,p) [args.dfdphan(x,p(1:end-3,:)) zeros(size(x,1),3,numel(t))];
    else
        data.dfdphan = @(t,x,p) [args.dfdphan(t,x,p(1:end-3,:)) zeros(size(x,1),3,numel(t))];
    end
end
if isempty(args.dfdthan)
    data.dfdthan = [];
else
    if data.autonomous
        data.dfdthan = @(t,x,p) args.dfdthan(x,p(1:end-3,:));
    else
        data.dfdthan = @(t,x,p) args.dfdthan(t,x,p(1:end-3,:));
    end
end
if ~data.autonomous && isempty(data.Om2idx)
    flag = any(strcmp(args.pnames,'Om2'));
    assert(flag, 'Om2 should include as a system parameter if its index is not specified');
    Om2idx = find(strcmp(args.pnames(:),'Om2'));
    data.Om2idx = Om2idx;
end
data.pnames  = args.pnames;

end