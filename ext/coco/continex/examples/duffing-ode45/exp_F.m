function [data y] = exp_F(data, x,p)

if isempty(data) || ~isfield(data, 'init')
    defaults.init          = true;
    defaults.sample_points = 50;
    defaults.TransTime     = 20;
    defaults.x0            = zeros(3,1);
    
    defaults.cdata.scale = 0.05;
    defaults.cdata.k1    = 0;
    defaults.cdata.k2    = 0;
    
    defaults.fcount      = 0;
    defaults.fcountMX    = 20;
    defaults.dev         = 0;
    
    data = coco_merge(defaults, data);

    TSample      = linspace(0,2*pi,data.sample_points+1);
    TSample      = TSample(1:end-1);
    data.Phi     = FTrans(TSample, data.sample_points);

    TSample      = TSample/(2*pi);
    data.TSample = TSample;
end

data.fcount = data.fcount + 1;
if data.fcount>=data.fcountMX
  data.fcount = 0;
  data.dev    = 0; %abs(1.0e-4*randn);
end

opts  = [];
opts  = odeset(opts, 'RelTol', 1.0e-5);
opts  = odeset(opts, 'Vectorized', 'on');
% if p(2)<0.13
%   opts  = odeset(opts, 'outputfcn', @odeplot);
% end

cdata    = data.cdata;
cdata.yr = @(t) IFTrans(x',p(2)*t);

assert(p(2)>=0.1, '%s: parameter ''om'' too small, aborting', mfilename);

x0      = data.x0;
T       = 2*pi/p(2);
[t,z]   = ode15s(@duff, [0 data.TransTime*T], x0, opts, p, cdata); %#ok<ASGLU>
x0      = z(end,:)';
[t,z]   = ode15s(@duff, [data.TSample 1]*T, x0, opts, p, cdata); %#ok<ASGLU>

data.x0 = z(end,:)';
y       = data.Phi*z(1:end-1,1);
dev     = data.dev;
y       = y.*(1 + dev*randn(size(y))) + dev*randn(size(y));
y       = y - x;
end

function y = IFTrans(c,t)
y = c*[
  ones(1,length(t))
  sin(t)
  sin(2*t)
  sin(3*t)
  cos(t)
  cos(2*t)
  cos(3*t)
  ];
end

function Phi = FTrans(t,N)
Phi = (2/N)*[
  ones(1,length(t))/2
  sin(t)
  sin(2*t)
  sin(3*t)
  cos(t)
  cos(2*t)
  cos(3*t)
  ];
end
