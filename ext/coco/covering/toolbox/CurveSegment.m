classdef CurveSegment < CurveSegmentBase & EventLocator
  
  properties (Access = private)
    % properties for interpolation polynomials
    A, B, C, D
    % function handle of interpolation function
    interpolate
    % class settings
    cont
  end
  
  properties (Access = public)
    % indicator for regular/singular initial point (default=regular)
    correct = true
    % projection condition
    prcond = struct('x', [], 'TS', [], 's', [], 'h', [])
  end
  
  methods (Access=private)
    
    function cseg = CurveSegment(opts, chart, init_flag)
      chart.rmProps = 'TS';
      cseg = cseg@CurveSegmentBase(opts, chart, init_flag);
      cseg.cont = cseg.class_props();
      switch cseg.cont.interp
        case 'cubic'
          cseg.interpolate = @CurveSegment.h3_interpolate;
        case 'linear'
          cseg.interpolate = @CurveSegment.lin_interpolate;
      end
    end
    
  end
  
  methods (Static)
    
    function [opts cseg] = create(opts, chart, prcond, x1)
      cseg              = CurveSegment(opts, chart, false);
      cseg.prcond       = prcond;
      cseg.ptlist       = { chart };
      cseg.curr_chart.x = x1;
      opts              = coco_emit(opts, 'update', cseg);
      % bug: enable optional recomputation of tangent vector
      % if cseg.correct
      %   [opts chart]    = cseg.update_TS(opts, chart);
      %   [opts chart]    = cseg.update_t(opts, chart);
      %   cseg.ptlist     = { chart };
      % end
      % opts              = coco_emit(opts, 'update2', cseg);
      % check residuum
      % [opts chart F]  = opts.efunc.F(opts, chart, chart.x);
    end
    
    function [opts cseg] = create_initial(opts, chart, varargin)
      cseg             = CurveSegment(opts, chart, true);
      cseg.interpolate = @CurveSegment.lin_interpolate;
      [opts cseg TS]   = cseg.init_TS(opts, varargin{:});
      cseg.prcond.x    = chart.x;
      cseg.prcond.TS   = TS;
      % bug: this should work, but leads to error in hermite interpolation
      % cseg.prcond.s    = zeros( size(TS,2), 1 );
      cseg.prcond.s    = [ 1 ; zeros( size(TS,2)-1, 1 ) ];
      cseg.prcond.h    = 0;
      opts             = coco_emit(opts, 'update', cseg);
    end
    
    function cont = class_props(S)
      persistent settings
      if isempty(settings) || (nargin>0 && ischar(S) && strcmpi(S, 'defaults'))
        settings = struct();
        settings.initTS = 'auto';
        settings.interp = 'cubic';
        settings.BP     = true;
      end
      if nargin>0 && isstruct(S)
        settings = coco_merge(settings, S);
      end
      cont = settings;
    end
    
    % function [opts data] = add_prcond(opts, dim)
    %   [opts data] = CurveSegmentBase.add_prcond(opts, dim);
    %   opts = coco_add_signal(opts, 'update2', mfilename);
    % end
  end
  
  methods (Static) % utility functions
    
    function [opts chart TS] = nullspace(opts, chart)
      % compute Jacobian of extended system
      fidx = coco_get_func_data(opts, 'cseg.prcond', 'fidx');
      [opts chart J] = opts.efunc.DFDX(opts, chart, chart.x);
      J(fidx,:) = [];
      
      % compute null space of linearisation of continuation problem
      % bug: move this to lsol and use lu as there (faster)
      [L U P]  = lu(J'); %#ok<ASGLU>
      [m n] = size(J);
      Y  = L(1:m, 1:m)' \ L(m+1:end, 1:m)';
      TS = P'*[Y; -speye(n-m)];
      TS = orth(full(TS));
    end
    
    function log_angle(opts, prio, LogLevel, t, idx, Name, StepName)
      for i=idx
        t2    = zeros(size(t));
        t2(i) = 1;
        coco_log(opts, prio, LogLevel, '%s: %s: angle(t,t%d) = % .2e[DEG]\n', ...
          Name, StepName, i, 180*subspace(t,t2)/pi);
      end
    end
    
    function Y = orth(X)
      Y = orth(X);
      % preserve handedness of base, this is necessary for guaranteeing
      % continuity of monitor functions that depend on the tangent space
      A = X' * Y; %#ok<PROP>
      Y(:,1) = sign(det(A))*Y(:,1); %#ok<PROP>
    end
    
    function chart = h3_interpolate(cseg, chart, la)
      [chart.x chart.t] = cseg.h3_eval(la);
    end
    
    function chart = lin_interpolate(cseg, chart, la)
      x0      = cseg.ptlist{  1}.x;
      x1      = cseg.ptlist{end}.x;
      t0      = cseg.ptlist{  1}.t;
      t1      = cseg.ptlist{end}.t;
      chart.x = (1-la)*x0 + la*x1;
      chart.t = (1-la)*t0 + la*t1; % bug: use interpolation on circle (?)
      chart.t = chart.t/norm(chart.t);
    end
    
    function [opts chart] = update_det(opts, chart)
      % compute Jacobian of extended system
      opts            = coco_emit(opts, 'set_mode', 'prcond');
      [opts chart J]  = opts.efunc.DFDX(opts, chart, chart.x);
      % % compute tangent space as solution of DF/DS*TS = [0...0 I]
      % fidx = coco_get_func_data(opts, 'cseg.prcond', 'fidx');
      % dim  = numel(fidx);
      % b    = zeros(size(chart.x,1), dim);
      % b(fidx,:)       = eye(dim,dim);
      [opts chart] = opts.lsol.solve(opts, chart, J, []);
    end
    
  end
  
  methods (Access=public) % interpolation functions
    
    function [opts, cseg, TS] = init_TS(cseg, opts, Valpha, h0, N)
      % compute initial tangent space matrix and, if ~cseg.correct,
      % cseg.curr_chart.TS.
      
      if nargin<3
        Valpha = 80;
      end
      if nargin<4
        h0 = 0.01;
      end
      if nargin<5
        N = 0;
      end
      
      chart = cseg.src_chart;
      dim   = numel(coco_get_func_data(opts, 'cseg.prcond', 'fidx'));
      if dim==0
        TS = zeros(size(chart.t,1),0);
        return;
      end
      
      idx = opts.efunc.p_idx;
      x0  = chart.x;
      t0  = chart.t;
      
      if all(t0 == 0)
        [opts, cseg.curr_chart, TS] = cseg.nullspace(opts, cseg.curr_chart);
        
        % find suitable vectors spanning parameter plane
        npars    = min(dim, numel(idx));
        dim_corr = max(1,sqrt(size(TS,1))-dim);
        t2       = [];
        for k=1:npars
          t1         = TS(:,k);
          % correct for dimensional scaling
          t1(idx(k)) = t1(idx(k))*dim_corr;
          t1         = t1(idx(k))/norm(t1);
          if Valpha>=180*acos(abs(t1))/pi
            t         = zeros(size(TS,1),1);
            t(idx(k)) = 1;
            t2 = [t2 t]; %#ok<AGROW>
          end
        end
        
        % complete new basis with "most orthogonal" elements from TS
        if ~isempty(t2)
          for i=(size(t2,2)+1):dim
            al = zeros(1,size(TS,2));
            for k=1:size(TS,2)
              al(k) = subspace(t2, TS(:,k));
            end
            [al, k] = max(al); %#ok<ASGLU>
            x       = TS(:,k)-t2*(t2'*TS(:,k));
            x       = x/norm(x);
            t2      = [t2 x]; %#ok<AGROW>
            TS(:,k) = [];
          end
          TS = t2;
        end
      else
        % initial direction provided (via coco_add_func), compute tangent
        % space using N+1 iterations of linearly implicit trapezoidal rule
        % N>=1 implies N correction steps on TS
        It  = 0;
        rep = true;
        while (rep && It<=N+1) || It<N+1
          It                = It+1;
          t0                = t0/norm(t0);
          chart.x           = x0 + h0*t0;
          [opts, chart, TS] = cseg.nullspace(opts, chart);
          if 180*subspace(TS,t0)/pi <= Valpha
            t1     = TS*(TS'*t0);
            t0     = t0 + t1/norm(t1);
            rep    = false;
          else
            [v, i] = max(abs(t0'*TS)); %#ok<ASGLU>
            t0     = t0 + TS(:,i); % tie break if TS (almost) orthogonal to t0
            rep    = true;         % perform at least one non-trivial iteration
          end
        end
        
        % make initial orientation of base vectors independent
        % of algorithm in CurveSegmentBase.nullspace
        % bug: in BPhandler
        % set sol.v such that parameter direction is preserved
%         for k=1:dim
%           TS(:,k) = sign(TS(:,k)'*t0)*TS(:,k);
%         end
        
        cseg.correct = false;
      end
      
      % make initial orientation of base vectors independent
      % of algorithm in CurveSegmentBase.nullspace
      for k=1:dim
        i = find(abs(TS(idx,k))>0, 1);
        if isempty(i)
          % go in direction of increasing norm
          x1 = x0+h0*TS(:,k);
          if norm(x1)<norm(x0)
            TS(:,k) = -TS(:,k);
          end
        else
          % go in direction of parameter
          TS(:,k) = sign(TS(idx(i),k))*TS(:,k);
        end
      end
      
      if ~cseg.correct
        cseg.curr_chart.TS = TS;
      end
      
    end
    
    function [opts chart] = update_TS(cseg, opts, chart)
      % compute Jacobian of extended system
      opts            = coco_emit(opts, 'set_mode', 'prcond');
      [opts chart J]  = opts.efunc.DFDX(opts, chart, chart.x);
      % compute tangent space as solution of DF/DS*TS = [0...0 I]
      fidx = coco_get_func_data(opts, 'cseg.prcond', 'fidx');
      dim  = numel(fidx);
      b    = zeros(size(chart.x,1), dim);
      b(fidx,:)       = eye(dim,dim);
      [opts chart TS] = opts.lsol.solve(opts, chart, J, b);
      chart.TS        = cseg.orth(TS);
    end
    
    function [opts chart] = update_t(cseg, opts, chart)
      % compute tangent vector at curve segment
      s       = (cseg.prcond.TS' * chart.TS) \ cseg.prcond.s;
      chart.t = chart.TS * s ;
      chart.t = chart.t/norm(chart.t);
    end
    
    function [opts cseg] = add_chart(cseg, opts, chart)
      if cseg.correct
        [opts chart]       = cseg.update_TS(opts, chart);
      end
      [opts chart]         = cseg.update_t(opts, chart);
      [opts chart chart.p] = opts.efunc.monitor_F(opts, chart, chart.x, chart.t);
      cseg.ptlist          = [ cseg.ptlist { chart } ];
    end
    
    function [opts cseg chart h] = chart_at(cseg, opts, la, evidx, varargin)
      chart = cseg.new_chart(opts, varargin{:});
      chart = cseg.interpolate(cseg, chart, la);
      
      [opts chart chart.p] = opts.efunc.monitor_F(opts, chart, chart.x, chart.t, evidx);
      
      x0 = cseg.prcond.x;
      TS = cseg.prcond.TS;
      s  = cseg.prcond.s;
      ss = s'*s; % bug: ss==1 always?
      h  = ((TS*s)'*(chart.x-x0))/ss;
    end
    
    function [opts cseg] = eval_p(cseg, opts, evidx)
      chart = cseg.curr_chart;
      [opts chart chart.p] = opts.efunc.monitor_F(opts, chart, chart.x, chart.t, evidx);
      cseg.curr_chart = chart;
    end
    
    function [opts chart] = update_p(cseg, opts, chart)
      opts.cseg            = cseg; % allow predict to call update_p
      [opts chart chart.p] = opts.efunc.monitor_F(opts, chart, chart.x, chart.t);
    end
    
    function [opts cseg idx] = insert_chart(cseg, opts, chart)
      x0  = cseg.prcond.x;
      TS  = cseg.prcond.TS;
      s   = cseg.prcond.s;
      t   = TS*s;
      h0  = t'*(chart.x-x0);
      h1  = cellfun(@(c) t'*(c.x-x0), cseg.ptlist(2:end));
      idx = find(h0<h1, 1)+1;
      
      if isempty(idx)
        % catch continuation events that landed on the wrong side of EP
        if h0-h1(end)<opts.corr.TOL
          idx = numel(cseg.ptlist);
        else
          coco_warn(opts, 1, opts.cont.LogLevel, ...
            '%s: point outside curve segment\n', mfilename);
          idx = numel(cseg.ptlist);
        end
      end
      
      if cseg.cont.BP
        [opts chart]       = cseg.update_det(opts, chart);
      end
      chart.TS             = TS;
      [opts chart chart.p] = opts.efunc.monitor_F(opts, chart, chart.x, chart.t);
      cseg.ptlist          = [ cseg.ptlist(1:idx-1) chart cseg.ptlist(idx:end) ];
    end
    
    function [opts cseg] = ev_init(cseg, opts)
      [opts cseg] = ev_init@EventLocator(cseg, opts);
      switch cseg.cont.interp
        case 'cubic'
          x0   = cseg.ptlist{  1}.x;
          x1   = cseg.ptlist{end}.x;
          t0   = cseg.ptlist{  1}.t;
          t1   = cseg.ptlist{end}.t;
          cseg = cseg.h3_icoeffs(x0, t0, x1, t1);
        case 'linear'
      end
    end
    
    function [opts cseg] = add_BP(cseg, opts, chart)
      [opts cseg] = add_BP@EventLocator(cseg, opts, chart);
      switch cseg.cont.interp
        case 'cubic'
          cseg = cseg.h3_adjust_icoeffs(cseg.u1, opts.cont.TOL);
        case 'linear'
      end
    end
    
  end
  
  methods (Access=private)
    
    function cseg = h3_icoeffs(cseg, x0, t0, x1, t1)
      %POLY_ICOEFFS Scale free Hermite interpolation coefficients.
      
      dx  = x1-x0;
      t00 = t0'*t0;
      t01 = t0'*t1;
      t11 = t1'*t1;
      dx0 = dx'*t0;
      dx1 = dx'*t1;
      
      % qubic spline with norm minimal nonlinear part
      % h(s) = C*s^2+D*s^3
      % Dk   = (d/ds)^k
      % sum_{i=0,1,2,3} ||(Dk h)(0)||^2 + ||(Dk h)(1)||^2 = Min!
      % the coefficients C and D depend on two scaling factors
      AA = [ 26*t00 23*t01 ; 23*t01 28*t11 ];
      bb = [ 49*dx0 ; 51*dx1 ];
      al = AA\bb;
      
      cseg.A = x0;
      cseg.B = al(1)*t0;
      cseg.D = al(2)*t1 + cseg.B - 2*dx;
      cseg.C = dx - cseg.B - cseg.D;
    end
    
    function cseg = h3_adjust_icoeffs(cseg, y, TOL)
      Rs  = norm(y-cseg.A);
      F0  = Rs;
      F1  = Rs-norm(cseg.h3_eval(1)-cseg.A);
      if abs(F1)<=TOL
        return
      end
      s0  = 0;
      s1  = 1;
      
      % compute first subdivision point
      if abs(F0)<=abs(F1)
        s = cseg.la2*s0+cseg.la1*s1;
      else
        s = cseg.la1*s0+cseg.la2*s1;
      end
      
      % locate point to remap
      while norm(s1-s0) >= TOL
        F = Rs-norm(cseg.h3_eval(s)-cseg.A);
        if F==0; break; end
        if F0*F<0
          s1 = s;
          F1 = F;
        else
          s0 = s;
          F0 = F;
        end
        if abs(F0)<=abs(F1)
          s = cseg.la2*s0+cseg.la1*s1;
        else
          s = cseg.la1*s0+cseg.la2*s1;
        end
      end
      
      % adjust B (rotate polynomial)
      cseg.B = (y-cseg.A)/s - s*(cseg.C + s*cseg.D);
      
      % rescale x-interval from [0 s] to [0 1]
      cseg.B = cseg.B*s;
      cseg.C = cseg.C*(s*s);
      cseg.D = cseg.D*(s*s*s);
    end
    
    function cseg = h3_xcoeffs(cseg, x0, t0, x1, t1)
      %POLY_XCOEFFS Scale free Hermite extrapolation coefficients.
      
      dx  = x1-x0;
      t00 = t0'*t0;
      t01 = t0'*t1;
      t11 = t1'*t1;
      dx0 = dx'*t0;
      dx1 = dx'*t1;
      
      % qubic spline with norm minimal nonlinear part
      % h(s) = C*s^2+D*s^3
      % Dk   = (d/ds)^k
      % sum_{i=0,1,2,3} ||(Dk h)(0)||^2 + ||(Dk h)(1)||^2 = Min!
      % the coefficients C and D depend on two scaling factors
      AA = [ 26*t00 23*t01 ; 23*t01 28*t11 ];
      bb = [ 49*dx0 ; 51*dx1 ];
      al = AA\bb;
      
      cseg.A = x1;
      cseg.B = al(2)*t1;
      cseg.D = al(1)*t0 + cseg.B - 2*dx;
      cseg.C = cseg.B + cseg.D - dx;
    end
    
    function [y t] = h3_eval(cseg, x)
      %H3_EVAL Evaluate Hermite polynomial and normalised tangent at X.
      
      rows = size(cseg.A,1);
      cols = numel(x);
      or   = ones(rows,1);
      oc   = ones(1,cols);
      
      x  = x(:)';
      x  = x(or,:);
      y  = cseg.A(:,oc) + x.*(cseg.B(:,oc) + x.*(cseg.C(:,oc) + x.*cseg.D(:,oc)));
      
      t  = cseg.B(:,oc) + x.*((2*cseg.C(:,oc)) + x.*(3*cseg.D(:,oc)));
      nt = sqrt(sum(t.*t,1));
      t  = t./nt(or,:);
    end
    
    function h3_plot(cseg, N, idx)
      
      if nargin<3
        [v i1] = sort( abs(cseg.B) ); %#ok<ASGLU>
        [v i2] = sort( max(abs(cseg.C),abs(cseg.D)) ); %#ok<ASGLU>
        if i1(end)~=i2(end)
          idx = [i1(end) i2(end)];
        else
          idx = [i1(end) i2(end-1)];
        end
      end
      
      x     = linspace(0, 1, N);
      [y t] = cseg.h3_eval(x);
      y     = y(idx,:);
      t     = t(idx,:);
      nt    = sqrt(sum(t.*t,1));
      or    = ones(size(t,1),1);
      t     = t./nt(or,:);
      hh    = 1.5*norm(y(:,1)-y(:,end))/N;
      
      switch numel(idx)
        
        case 2
          plot(y(1,:), y(2,:), 'b.-');
          hold on
          for i=1:N
            xx = [y(1,i) y(1,i)+hh*t(1,i)];
            yy = [y(2,i) y(2,i)+hh*t(2,i)];
            plot(xx, yy, 'r-');
          end
          hold off
          
        case 3
          plot3(y(1,:), y(2,:), y(3,:), 'b.-');
          hold on
          for i=1:N
            xx = [y(1,i) y(1,i)+hh*t(1,i)];
            yy = [y(2,i) y(2,i)+hh*t(2,i)];
            zz = [y(3,i) y(3,i)+hh*t(3,i)];
            plot3(xx, yy, zz, 'r-');
          end
          hold off
          
      end
      
    end
    
  end
  
  methods (Access = public)
    
    function [opts cseg idx] = add_SP(cseg, opts, chart)
      [opts chart]    = opts.cseg.update_det(opts, chart); % Harry added: Compute determinant of restricted problem
      [opts cseg idx] = add_SP@EventLocator(cseg, opts, chart);
    end
    
  end
  
end
