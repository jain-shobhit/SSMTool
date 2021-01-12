classdef cseg_continex < CurveSegmentBase & EventLocator
  % cseg_continex   Basic curve segment class with non-linear covering condition.
  
  properties (Access = public)
    prcond     = struct('x', [], 'TS', [], 's', [], 'h', [])
    correct    = true
    trial      = 1;
    corr       = struct();
    MaxAbsDist = [];
    base_chart = [];
  end
  
  methods (Access=private)
    
    function cseg = cseg_continex(prob, chart, init_flag)
      cseg = cseg@CurveSegmentBase(prob, chart, init_flag);
    end
    
  end
  
  methods (Static)
    
    function cont = get_settings(cont)
      defaults.shape = 'ellipsis';  % shape class for arc-length condition
      cont           = coco_merge(defaults, cont);
      if ischar(cont.shape)
        shclass    = sprintf('continex_shape_%s', cont.shape);
        cont.shape = str2func(shclass);
      else
        shclass    = func2str(cont.shape);
      end
      func = str2func(sprintf('%s.get_settings', shclass));
      cont = func(cont);
    end
    
    function [prob cseg] = create_initial(prob, chart, varargin)
      cseg            = cseg_continex(prob, chart, true);
      cseg.base_chart = chart;
      cseg.prcond.x   = chart.x;
      cseg.prcond.TS  = cseg.init_TS(prob);
      cseg.prcond.s   = 1;
      cseg.prcond.h   = 0;
      prob            = coco_emit(prob, 'update', cseg);
      prob            = coco_emit(prob, 'set_mode', 'prcond');
    end
    
    function [prob cseg] = create(prob, chart, prcond)
      cseg            = cseg_continex(prob, chart, false);
      cseg.base_chart = chart;
      cseg.prcond     = prcond;
      cseg.ptlist     = { chart };
      
      data = coco_get_func_data(prob, 'cseg.prcond', 'data');
      
      if chart.pt == 0
        mode = 'prcond';
        h    = prcond.h;
      else
        mode = 'sh_cond';
        h    = data.shape.h;
        cseg.MaxAbsDist = data.shape.scale*data.shape.MaxAbsDist;
      end
      
      x1 = chart.x + (data.shape.scale*h)*(prcond.TS*prcond.s);
      cseg.curr_chart.x = x1;
      
      prob = coco_emit(prob, 'update', cseg);
      prob = coco_emit(prob, 'set_mode', mode);
      
      % if chart.pt>0
      %   % [data f] = cseg.prcond_F(prob, data, x1);
      %   % coco_print(prob, 3, '%s: scale=%.2e, h=%.2e, prcond_F=% .2e\n', ...
      %   %   mfilename, data.shape.scale, data.shape.scale*prcond.h, f);
      %   coco_print(prob, 3, '%s: scale=%.2e, h=%.2e\n', ...
      %     mfilename, data.shape.scale, data.shape.scale*h);
      % end
      
    end
    
  end
  
  methods
    
    function [prob cseg] = reset(cseg, prob)
      cseg.src_chart  = cseg.base_chart;
      cseg.curr_chart = cseg.init_chart_from(prob, cseg.src_chart);
      
      data = coco_get_func_data(prob, 'cseg.prcond', 'data');
      
      if cseg.base_chart.pt == 0
        mode = 'prcond';
        h    = cseg.prcond.h;
      else
        mode = 'sh_cond';
        h    = data.shape.h;
        cseg.MaxAbsDist = data.shape.scale*data.shape.MaxAbsDist;
      end
      
      x1 = cseg.base_chart.x + (data.shape.scale*h)*(cseg.prcond.TS*cseg.prcond.s);
      cseg.curr_chart.x = x1;
      
      prob = coco_emit(prob, 'update', cseg);
      prob = coco_emit(prob, 'set_mode', mode);
      
      % % [data f] = cseg.prcond_F(prob, data, x1);
      % % coco_print(prob, 3, '%s: scale=%.2e, h=%.2e, prcond_F=% .2e\n', ...
      % %   mfilename, data.shape.scale, data.shape.scale*cseg.prcond.h, f);
      % coco_print(prob, 3, '%s: scale=%.2e, h=%.2e\n', ...
      %   mfilename, data.shape.scale, data.shape.scale*cseg.prcond.h);
      
    end
    
    function TS = init_TS(cseg, prob)
      pidx = coco_get_func_data(prob, 'efunc', 'pidx');
      assert(numel(pidx) >= 1, ...
        '%s: no continuation parameter present', mfilename);
      % go in direction of first active parameter
      TS          = zeros(size(cseg.src_chart.x));
      TS(pidx(1)) = 1;
    end
    
    function [prob chart] = update_TS_and_t(cseg, prob, chart, pts)
      if isempty(cseg.base_chart) || chart.pt==0
        % use tangent space from projection condition
        t  = cseg.prcond.TS*chart.s;
      elseif nargin<4 || size(pts,2)<2
        % compute normalised secant direction
        t = (chart.x - cseg.base_chart.x);
        t = t/norm(t);
      else
        % compute least squares fit of direction
        [U S]    = eig(cov(pts')');
        [ss idx] = sort(diag(S), 'descend'); %#ok<ASGLU>
        t        = U(:,idx(1));
        t2       = (chart.x - cseg.base_chart.x);
        t2       = t2/norm(t2);
        if t'*t2 < 0 % align with secant direction
          t = -t;
        end
      end
      chart.TS = chart.s*t;
      chart.t  = t;
    end
    
    function [prob cseg] = add_chart(cseg, prob, chart, varargin)
      [prob chart]         = cseg.update_TS_and_t(prob, chart, varargin{:});
      % Evaluate monitor functions at this point only, use
      % isfield(chart,'monitor') to check if evaluation occurs here.
      chart.monitor        = true;
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t);
      chart                = rmfield(chart, 'monitor');
      cseg.ptlist          = [ cseg.ptlist { chart } ];
    end
    
    function [prob cseg] = ev_init(cseg, prob)
      [prob cseg] = ev_init@EventLocator(cseg, prob);
      u0 = cseg.ptlist{1}.x;
      u1 = cseg.ptlist{end}.x;
      t  = (u1-u0);
      h  = norm(t);
      t  = t/h;
      cseg.prcond.x  = u0;
      cseg.prcond.TS = t*cseg.ptlist{1}.s;
      cseg.prcond.s  = cseg.ptlist{1}.s;
      cseg.prcond.h  = h;
      prob = coco_emit(prob, 'update', cseg);
    end
    
    function [prob cseg chart h] = chart_at(cseg, prob, la, evidx, varargin)
      chart = cseg.new_chart(prob, varargin{:});
      
      % bug: interpolate monitor functions?
      x0      = cseg.ptlist{  1}.x;
      x1      = cseg.ptlist{end}.x;
      % p0      = cseg.ptlist{  1}.p;
      % p1      = cseg.ptlist{end}.p;
      t0      = cseg.ptlist{  1}.t;
      t1      = cseg.ptlist{end}.t;
      chart.x = (1-la)*x0 + la*x1;
      % chart.p = (1-la)*p0 + la*p1;
      chart.t = (1-la)*t0 + la*t1; % bug: use interpolation on circle (?)
      chart.t = chart.t/norm(chart.t);
      
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t, evidx);
      
      x0 = cseg.prcond.x;
      TS = cseg.prcond.TS;
      s  = cseg.prcond.s;
      ss = s'*s; % bug: ss==1 always?
      h  = ((TS*s)'*(chart.x-x0))/ss;
    end

    function [prob cseg] = eval_p(cseg, prob, evidx)
      chart = cseg.curr_chart;
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t, evidx);
      cseg.curr_chart = chart;
    end
    
    function [prob cseg idx] = insert_chart(cseg, prob, chart)
      x0  = cseg.prcond.x;
      TS  = cseg.prcond.TS;
      s   = cseg.prcond.s;
      t   = TS*s;
      h0  = t'*(chart.x-x0);
      h1  = cellfun(@(c) t'*(c.x-x0), cseg.ptlist(2:end));
      idx = find(h0<h1, 1)+1;
      
      if isempty(idx)
        % catch continuation events that landed on the wrong side of EP
        if h0-h1(end)<prob.corr.TOL
          idx = numel(cseg.ptlist);
        else
          error('%s: point outside curve segment', mfilename);
        end
      end
      
      chart.TS             = TS;
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t);
      cseg.ptlist          = [ cseg.ptlist(1:idx-1) chart cseg.ptlist(idx:end) ];
    end
    
    function plot_prcond(cseg, prob, phan, xidx, yidx) %#ok<MANU>
      data = coco_get_func_data(prob, 'cseg.prcond', 'data');
      switch data.mode
        case 1
          data.shape.plot_shape(phan, xidx, yidx);
        case 3
          % plot tangent to shape in (x,||y||) coordinates
          % t = data.TS'*data.s;
          % u = data.u0;
          % h = cseg.prcond.h;
          % t = t/norm(t([xidx yidx]));
          % x = [u(xidx)-h*t(xidx) u(xidx)+h*t(xidx)];
          % y = [u(yidx)-h*t(yidx) u(yidx)+h*t(yidx)];
          % plot(phan, x,y,'y-', u(xidx),u(yidx),'bp');
      end
    end
    
    function plot_point(cseg, prob, phan, x, xidx, yidx, varargin) %#ok<MANU>
      data = coco_get_func_data(prob, 'cseg.prcond', 'data');
      switch data.mode
        case 1
          data.shape.plot_point(phan, x, xidx, yidx, varargin{:});
        case 3
          % plot points on line tangent to shape in (x,||y||) coordinates
          % t = data.TS'*data.s;
          % u = data.u0;
          % h = cseg.prcond.h;
          % t = t/norm(t([xidx yidx]));
          % x = [u(xidx)-h*t(xidx) u(xidx)+h*t(xidx)];
          % y = [u(yidx)-h*t(yidx) u(yidx)+h*t(yidx)];
          % plot(phan, x,y,'y-', u(xidx),u(yidx),'bp');
      end
    end
    
  end
  
  methods (Static)
    
    function prob = add_prcond(prob, cont, dim)
      prob = coco_add_signal(prob, 'update',    mfilename);
      prob = coco_add_signal(prob, 'set_mode',  mfilename);
      prob = coco_add_signal(prob, 'update_h',  mfilename);
      prob = coco_add_signal(prob, 'fix_mfunc', mfilename);
      prob = coco_add_signal(prob, 'shape_set_scale', mfilename);
      
      data = struct('dim', dim, 'mode', 0, ...
        'u0', [], 'TS', [], 's', [], 'h', 0);
      
      data.shape = cont.shape(cont);
      
      data = coco_func_data(data);
      prob = coco_add_func(prob, 'cseg.prcond', ...
        @cseg_continex.prcond_F, @cseg_continex.prcond_DFDX, ...
        data, 'zero', 'xidx', 'all', 'fdim', 1);
      
      prob = coco_add_slot(prob, 'cseg.update', ...
        @cseg_continex.update, data, 'update');
      prob = coco_add_slot(prob, 'cseg.set_mode', ...
        @cseg_continex.set_mode, data, 'set_mode');
      prob = coco_add_slot(prob, 'cseg.update_h', ...
        @cseg_continex.update_h, data, 'update_h');
      prob = coco_add_slot(prob, 'cseg.fix_mfunc', ...
        @cseg_continex.fix_mfunc, data, 'fix_mfunc');
      
      prob = coco_add_slot(prob, 'cseg.set_scale', ...
        @cseg_continex.set_scale, data, 'shape_set_scale');
      prob = coco_add_slot(prob, 'cseg.set_sample', ...
        @cseg_continex.set_sample, data, 'corr_sample');
      prob = coco_add_slot(prob, 'cseg.pull_back', ...
        @cseg_continex.pull_back, data, 'corr_pull_back');

      prob = coco_add_slot(prob, 'cseg.remesh', ...
        @cseg_continex.remesh, data, 'remesh');
    end
    
    function [data f] = prcond_F(prob, data, u)  %#ok<INUSD>
      
      switch data.mode
        
        case 0 % initialise f for coco_add_func
          f = zeros(data.dim,1);
          
        case 1 % non-linear shape condition
          f = data.shape.F(u);
          
        case 2 % fix embedded test function along curve segment
          f = u(data.fixpar_idx) - data.fixpar_val;
          TSjj = data.TS(data.jj,:);
          for j=data.j
            f = [f; (data.TS(j,:)-data.ss(j)*TSjj)*(u-data.u0)]; %#ok<AGROW>
          end
          
        case 3 % linear projection condition for event location
          f = data.TS*(u-data.u0) - data.s*data.h;
          
      end
      
    end
    
    function [data J] = prcond_DFDX(prob, data, u) %#ok<INUSD>
      
      switch data.mode
        
        case 0
          J = sparse(data.dim, numel(u));
          
        case 1
          J = data.shape.DFDU(u);
          % [data J2] = coco_ezDFDX('f(o,d,x)', prob, data, @cseg_continex.prcond_F, u);
          
        case 2
          J    = sparse(1, data.fixpar_idx, 1, 1, numel(u));
          TSjj = data.TS(data.jj,:);
          for j=data.j
            J = [ J ; data.TS(j,:)-data.ss(j)*TSjj ]; %#ok<AGROW>
          end
          
        case 3
          J = data.TS;
          
      end
      
    end
    
    function data = update(prob, data, cseg, varargin) %#ok<INUSD>
      data.u0 = cseg.prcond.x;
      data.TS = cseg.prcond.TS';
      data.s  = cseg.prcond.s;
      data.h  = cseg.prcond.h;

      [data.ss data.jj] = max(abs(data.s));
      data.ss           = data.s/data.s(data.jj);
      data.j            = setdiff(1:data.dim, data.jj);
      
      data.shape = data.shape.update(cseg.prcond);
      
      data.mode = 3;
    end
    
    function data = set_mode(prob, data, mode, varargin) %#ok<INUSD>
      switch lower(mode)
        case {'init' 'check_res'}
          data.mode = 0;
        case 'sh_cond'
          data.mode = 1;
        case 'fix_mfunc'
          data.mode = 2;
        case 'prcond'
          data.mode = 3;
      end
    end
    
    function data = update_h(prob, data, h, varargin) %#ok<INUSD>
      data.h    = h;
      data.mode = 3;
    end
    
    function data = fix_mfunc(prob, data, mf_idx, mf_val, varargin) %#ok<INUSD>
      data.fixpar_idx = mf_idx;
      data.fixpar_val = mf_val;
      data.mode       = 2;
    end
    
    function data = set_scale(prob, data, scale, varargin) %#ok<INUSD>
      data.shape = data.shape.set_scale(scale);
    end
    
    function data = set_sample(prob, data, u, varargin)
      if data.mode==1
        u0        = data.shape.pull_back(u);
        J         = data.shape.DFDU(u0);
        t         = data.shape.t;
        ts        = (u0-data.shape.u1);
        ts        = ts/norm(ts);
        data.u0   = u0;
        data.TS   = J/norm(J);
        data.s    = 1;
        data.h    = 0;
        data.mode = 3;
        coco_log(prob, 1, prob.cont.LogLevel, '%s: angle(t ,ts) = %7.2f\n', ...
          mfilename, real(acos(t'*ts))*180/pi);
        coco_log(prob, 1, prob.cont.LogLevel, '%s: angle(ts,TS) = %7.2f\n', ...
          mfilename, real(acos(ts'*data.TS'))*180/pi);
        coco_log(prob, 1, prob.cont.LogLevel, '%s: angle(t ,TS) = %7.2f\n', ...
          mfilename, real(acos(t'*data.TS'))*180/pi);
      end
    end
    
    function [data u] = pull_back(prob, data, u, varargin) %#ok<INUSD>
      if data.mode==1
        u = data.shape.pull_back(u);
      end
    end
    
    function data = remesh(prob, data, varargin) %#ok<INUSD>
      % set projection condition back to init mode so that change_func
      % succeeds when called in coco_remesh
      data.mode = 0;
    end
    
  end
  
end
