classdef continex_shape_parabola < continex_shape
  
  properties (Access=public)
  end
  
  methods (Access=public)
    
    function shape = continex_shape_parabola(cont)
      shape = shape@continex_shape(cont);
      
      % construct parabola h*(x-h)+y^2=0
      shape.h = cont.h_max;
      
      shape.MaxAbsDist = sqrt(6)*shape.h;
    end
    
    function f = F(shape, u)
      x = shape.t'*(u-shape.u1);
      y = (u - shape.u1) - shape.t*x;
      f = shape.scale*shape.h*(x-shape.scale*shape.h) + y'*y;
    end
    
    function J = DFDU(shape, u)
      x = shape.t'*(u-shape.u1);
      y = (u - shape.u1) - shape.t*x;
      J = shape.scale*shape.h*shape.t' + 2*(y'-(y'*shape.t)*shape.t');
    end
    
    function shape = update(shape, prcond)
      shape.t  = prcond.TS*prcond.s;
      shape.u1 = prcond.x;
    end
    
    function x = pull_back(shape, x)
      ts = x-shape.u1;
      xx = shape.t'*ts;
      yy = ts - shape.t*xx;
      sh = shape.scale*shape.h;
      a  = 0.5*xx/sh;
      b  = (yy'*yy)/sh^2;
      r  = 1/(a+sqrt(a^2+b));
      x  = shape.u1+r*ts;
    end
    
  end
  
  methods (Access=public,Static=true)
    
    function cont = get_settings(cont)
      cont = continex_shape.get_settings(cont);
    end
    
  end
  
  methods (Access=public)
    
    function plot_shape(shape, phan, xidx, yidx)
      sh = shape.scale*shape.h;
      sr = shape.scale*shape.MaxAbsDist;
      y = linspace(-sqrt(2)*sh, sqrt(2)*sh, 100);
      T = shape.t([xidx yidx]);
      T = T/norm(T);
      N = [T(2) ; -T(1)];
      x = kron(sh-y.^2/sh, T) + kron(y, N);
      u = x + repmat(shape.u1([xidx yidx]), size(y));
      plot(phan, shape.u1(xidx), shape.u1(yidx), 'b*');
      plot(phan, u(1,:), u(2,:), 'b-');
      s = linspace(12*pi/16, 20*pi/16, 30);
      x = kron(sh+sr*cos(s), T) + kron(sr*sin(s), N);
      u = x + repmat(shape.u1([xidx yidx]), size(s));
      plot(phan, u(1,:), u(2,:), 'c-');
    end
    
    function plot_point(shape, phan, u, xidx, yidx, sample) %#ok<INUSD>
      T = shape.t([xidx yidx]);
      T = T/norm(T);
      N = [T(2) ; -T(1)];
      x = shape.t'*(u-shape.u1);
      y = (u - shape.u1) - shape.t*x;
      s = sign(y([xidx yidx])'*N);
      U = shape.u1([xidx yidx])+x*T+(s*norm(y))*N;
      plot(phan, U(1), U(2), 'mo');
    end
    
  end
    
end
