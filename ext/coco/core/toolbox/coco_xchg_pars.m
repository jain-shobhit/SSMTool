function opts = coco_xchg_pars(varargin)
%COCO_XCHG_PARS   Exchange inactive, active and internal parameters.
%
%   OPTS = COCO_XCHG_PARS([OPTS], PAR1, PAR2) exchanges parameter PAR1 with
%   parameter PAR2. PAR1 and PAR2 are strings containing the names of the
%   parameters to exchange. After the exchange, PAR1 will play the role of
%   PAR2 and vice versa in a subsequent continuation.
%
%   Note that internal parameters can be exchanged automatically by
%   parameter overspecification.

%% affected fields in opts
%
%    opts.xchg - cell array containing a list of pairs of parameter
%                names to exchange in coco_close_efunc

%% check for input argument opts
%  varargin = { [opts], PAR1, PAR2 }

argidx = 1;
if isempty(varargin{argidx}) || isstruct(varargin{argidx})
	opts   = varargin{argidx};
	argidx = argidx + 1;
else
	opts = [];
end

%% parse input arguments

par1 = varargin{argidx  };
par2 = varargin{argidx+1};

%% add entry to exchange list

if ~isfield(opts, 'efunc')
  opts.efunc = [];
end

if isfield(opts.efunc, 'xchg')
  opts.efunc.xchg = [opts.efunc.xchg ; {par1 par2}];
else
	opts.efunc.xchg = {par1 par2};
end

