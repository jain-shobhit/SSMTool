classdef continex_shape_ellipsis < continex_shape
  
  properties (Access=public)
    a,b,e, A,B,C
  end
  
  methods (Access=public)
    
    function shape = continex_shape_ellipsis(cont)
      shape = shape@continex_shape(cont);
      
      % construct ellipsoid with e = |ela*a| and (1+ela)*a = h
      shape.h = cont.h_max;
      shape.a = shape.h/(1+cont.ela);
      shape.b = shape.a*sqrt(1-cont.ela^2);
      shape.e = cont.ela*shape.a;
      
      shape.A = shape.b/shape.a;
      shape.B = shape.a/shape.b;
      shape.C = shape.a*shape.b;
      
      shape.MaxAbsDist = shape.a+abs(shape.e);
    end
    
    function f = F(shape, u)
      x = shape.t*(shape.t'*(u-shape.u1));
      y = (u - shape.u1) - x;
      f = sqrt(shape.A*(x'*x) + shape.B*(y'*y)) - shape.scale*sqrt(shape.C);
    end
    
    function J = DFDU(shape, u)
      x = shape.t*(shape.t'*(u-shape.u1));
      y = (u - shape.u1) - x;
      J = shape.A*((x'*shape.t)*shape.t') + shape.B*(y'-(y'*shape.t)*shape.t');
      J = J/sqrt(shape.A*(x'*x) + shape.B*(y'*y));
    end
    
    function shape = update(shape, prcond)
      shape.t  = prcond.TS*prcond.s;
      shape.u1 = prcond.x + shape.scale*shape.e*shape.t;
    end
    
    function x = pull_back(shape, x)
      ts = x-shape.u1;
      ts = ts/norm(ts);
      D2 = (shape.t'*ts)^2;
      r  = sqrt(shape.C)/sqrt(shape.A*D2+shape.B*(1-D2));
      x  = shape.u1+shape.scale*r*ts;
    end
    
  end
  
  methods (Access=public,Static=true)
    
    function cont = get_settings(cont)
      cont         = continex_shape.get_settings(cont);
      gmi          = 2/(1+sqrt(5)); % inverse golden mean
      defaults.ela = gmi; % excentricity factor e = |ela*a|
      cont         = coco_merge(defaults, cont);
      if ischar(cont.ela)
        switch lower(cont.ela)
          case { '+', '+g', 'g' }
            cont.ela = gmi;
          case { '-', '-g' }
            cont.ela = -gmi;
          otherwise
            error('%s: unknown setting ''%s'' for property ''ela''', ...
              mfilename, cont.ela);
        end
      end
    end
    
  end
  
  methods (Access=public)
    
    function plot_shape(shape, phan, xidx, yidx)
      sr = shape.scale*shape.MaxAbsDist;
      s = linspace(0, 2*pi, 100);
      T = shape.t([xidx yidx]);
      T = T/norm(T);
      N = [T(2) ; -T(1)];
      x = kron(shape.scale*shape.a*cos(s), T) + kron(shape.scale*shape.b*sin(s), N);
      u = x + repmat(shape.u1([xidx yidx]), size(s));
      plot(phan, shape.u1(xidx), shape.u1(yidx), 'b*');
      plot(phan, shape.u1(xidx)-shape.scale*shape.e*T(1), shape.u1(yidx)-shape.scale*shape.e*T(2), 'b*');
      plot(phan, u(1,:), u(2,:), 'b-');
      s = linspace(13*pi/16, 19*pi/16, 30);
      x = kron(shape.scale*shape.a+sr*cos(s), T) + kron(sr*sin(s), N);
      u = x + repmat(shape.u1([xidx yidx]), size(s));
      plot(phan, u(1,:), u(2,:), 'c-');
    end
    
    function plot_point(shape, phan, u, xidx, yidx, sample) %#ok<INUSD>
      T = shape.t([xidx yidx]);
      T = T/norm(T);
      N = [T(2) ; -T(1)];
      x = shape.t*(shape.t'*(u-shape.u1));
      y = (u - shape.u1) - x;
      s = sign(y([xidx yidx])'*N);
      U = shape.u1([xidx yidx])+(shape.t'*x)*T+(s*norm(y))*N;
      plot(phan, U(1), U(2), 'mo');
    end
    
  end
    
end
