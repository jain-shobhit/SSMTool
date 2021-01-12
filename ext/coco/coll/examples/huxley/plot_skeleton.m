function plot_skeleton(p0, vs, vu, varargin)
vs = sign(vs(1))*vs;
vu = sign(vu(1))*vu;

opts = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-6, 'NormControl', 'on');

X = linspace(-0.3,0,300);
[t, x] = ode45(@(t,x) huxley(x,p0), [0 11.5], -1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t, x] = ode45(@(t,x) -huxley(x,p0), [0 11.5], -1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

X = linspace(0,0.5,500);
[t, x] = ode45(@(t,x) huxley(x,p0), [0 13.5], 1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t, x] = ode45(@(t,x) -huxley(x,p0), [0 13.5], 1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([0;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

X = linspace(1,1.3,300);
[t, x] = ode45(@(t,x) huxley(x,p0), [0 11.5], [1;0]+1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t, x] = ode45(@(t,x) -huxley(x,p0), [0 11.5], [1;0]+1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

X = linspace(0.5,1,500);
[t, x] = ode45(@(t,x) huxley(x,p0), [0 13.5], [1;0]-1.0e-4*vu, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})
[t, x] = ode45(@(t,x) -huxley(x,p0), [0 13.5], [1;0]-1.0e-4*vs, opts); %#ok<ASGLU>
Y = interp1([1;x(:,1)], [0;x(:,2)], X);
plot(X,Y,varargin{:})

end
