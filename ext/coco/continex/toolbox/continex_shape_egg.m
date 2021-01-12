classdef continex_shape_egg < continex_shape
  
  properties (Access=public)
    % c = [ 14 -24 6 4 ]/27
    r  = @(x) 1 + (x<1).*( (14/27)-(24/27)*x+( 6/27)*x.^2+(4/27)*x.^3);
    rs = @(x)     (x<1).*(-(24/27)+(12/27)*x+(12/27)*x.^2);
  end
  
  methods (Access=public)
    
    function shape = continex_shape_egg(cont)
      shape = shape@continex_shape(cont);
      
      % construct eggshape sqrt(x^2+y^2) = h*r(x)
      shape.h = cont.h_max;
      
      shape.MaxAbsDist = 3.32*shape.h;
    end
    
    function f = F(shape, u)
      x  = shape.t'*(u-shape.u1);
      y  = (u - shape.u1) - shape.t*x;
      z  = u-shape.u1;
      sh = shape.scale*shape.h;
      f  = sqrt(z'*z) - sh*shape.r(x/sh);
    end
    
    function J = DFDU(shape, u)
      x  = shape.t'*(u-shape.u1);
      y  = (u - shape.u1) - shape.t*x;
      z  = u-shape.u1;
      sh = shape.scale*shape.h;
      J  = z'/sqrt(z'*z) - shape.rs(x/sh)*shape.t';
    end
    
    function shape = update(shape, prcond)
      shape.t  = prcond.TS*prcond.s;
      shape.u1 = prcond.x;
    end
    
    function x = pull_back(shape, x)
      r1 = 0.9;
      r2 = 3.1;
      sh = shape.scale*shape.h;
      ts = x-shape.u1;
      rr = norm(ts)/sh;
      ts = ts/rr;
      g1 = 2/(1+sqrt(5)); % 0.6180...
      g2 = 1-g1;          % 0.3820...
      x1 = shape.u1+r1*ts;
      x2 = shape.u1+r2*ts;
      f1 = shape.F(x1);
      f2 = shape.F(x2);
      ff = shape.F(x);
      if r1<rr && rr<r2
        if f1*ff<=0
          r2 = rr;
          f2 = ff;
        else
          r1 = rr;
          f1 = ff;
        end
      end
      % fprintf('pull_back: (F1,F,F2) = (% .2e,% .2e,% .2e)\n', f1,ff,f2);
      while r2-r1>1.0e-4
        if abs(f1)<=abs(f2)
          rr = g1*r1+g2*r2;
        else
          rr = g2*r1+g1*r2;
        end
        x  = shape.u1+rr*ts;
        ff = shape.F(x);
        if f1*ff<=0
          r2 = rr;
          f2 = ff;
        else
          r1 = rr;
          f1 = ff;
        end
        % fprintf('pull_back: (F1,F,F2) = (% .2e,% .2e,% .2e)\n', f1,ff,f2);
      end
      % fprintf('pull_back: done\n');
    end
    
    function x = pull_back_old(shape, x)
      % fprintf('pull_back: F = % .2e\n', shape.F(x));
      while abs(shape.F(x))>1.0e-10
        ts = x-shape.u1;
        xx = shape.t'*ts;
        yy = ts - shape.t*xx;
        sh = shape.scale*shape.h;
        rr = sh*shape.r(xx/sh)/norm(ts);
        x  = shape.u1+rr*ts;
        % fprintf('pull_back: F = % .2e\n', shape.F(x));
      end
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
      x  = linspace(-2*sh,sh, 100);
      y  = sqrt((sh*shape.r(x/sh)).^2-x.^2);
      T = shape.t([xidx yidx]);
      T = T/norm(T);
      N = [T(2) ; -T(1)];
      plot(phan, shape.u1(xidx), shape.u1(yidx), 'b*');
      u = kron(x, T) + kron(y, N);
      u = u + repmat(shape.u1([xidx yidx]), size(y));
      plot(phan, u(1,:), u(2,:), 'b-', u(1,:), u(2,:), 'b-');
      u = kron(x, T) - kron(y, N);
      u = u + repmat(shape.u1([xidx yidx]), size(y));
      plot(phan, u(1,:), u(2,:), 'b-', u(1,:), u(2,:), 'b-');
      s = linspace(11*pi/16, 21*pi/16, 30);
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
