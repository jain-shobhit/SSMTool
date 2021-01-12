classdef atlas2_x < AtlasBase
  
  properties (Access = private)
    cont     = struct()
    boundary = []
    charts   = {}
    next_pt  = -1;
    next_id  = 1;
  end
  
  methods ( Access = public )
    
    function atlas = atlas2_x(prob, cont, dim)
      assert(dim==2, '%s: wrong manifold dimension', mfilename);
      atlas      = atlas@AtlasBase(prob);
      atlas.cont = atlas.get_settings(cont);
    end
    
  end
  
  methods (Static)
    
    function [prob cont atlas] = create(prob, cont, dim)
      atlas = atlas2_x(prob, cont, dim);
      prob  = CurveSegment.add_prcond(prob, dim);
      prob  = coco_add_slot(prob, 'atlas', @atlas.save_atlas, [], 'save_bd');
    end
    
    function [data res] = save_atlas(opts, data, varargin)
%       if isempty(opts.atlas.boundary)
        res = opts.atlas.charts;
%       else
%         res = {};
%       end
    end
    
  end
  
  methods (Static, Access = private) % settings
    
    function cont = get_settings(cont)
      defaults.h     = 0.1 ; % continuation step size
      defaults.PtMX  = 100 ; % number of continuation steps
      defaults.theta = 0.5 ; % step size of theta method
      defaults.almax = 15  ; % maximum angle between tangent spaces
      defaults.Rmarg = 0.95; % boundary of margin of width=R*(1-Rmarg)
      defaults.Ndirs = 6   ; % number of initial directions
      cont           = coco_merge(defaults, cont);
      cont.PtMX      = abs(cont.PtMX);
      cont.almax     = cont.almax*pi/180;
      cont.Ndirs     = max(3, ceil(cont.Ndirs));
      al             = ((1:cont.Ndirs)-1)*(2*pi/cont.Ndirs);
      r1             = cont.h * sqrt(1.1+tan(pi/cont.Ndirs)^2);
      cont.s0        = [cos(al);sin(al)];     % initial directions
      cont.v0        = r1*ones(cont.Ndirs,1); % initial vertices
      cont.bv0       = 1:cont.Ndirs;          % list of boundary vertices
      cont.nb0       = zeros(1,cont.Ndirs);   % list of neighbours
    end
    
  end
  
  methods % interface methods
    
    function [prob atlas cseg correct] = init_prcond(atlas, prob, chart)
      [atlas chart] = atlas.set_fields(chart);
      chart.R       = 0;
      chart.bv      = [];
      chart.pt_type = 'IP';
      chart.ep_flag = 1;
      [prob cseg]   = CurveSegment.create_initial(prob, chart);
      correct       = cseg.correct;
    end
    
    function [prob atlas cseg] = init_chart(atlas, prob, cseg)
      [atlas chart]   = atlas.set_fields(cseg.curr_chart);
      [prob cseg]     = cseg.add_chart(prob, chart);
      cseg.curr_chart = cseg.ptlist{1};
    end
    
    function [opts atlas cseg] = init_admissible(atlas, opts, cseg, S)
      if S.ep_flag == 1
        % compute subset of admissible directions
        flags           = ~(S.dir_flags==atlas.IsAdmissible);
        chart           = cseg.curr_chart;
        chart.v(flags)  = 0.5*(atlas.cont.Rmarg*chart.R);
        chart.bv(flags) = [];
        cseg.curr_chart = chart;
        if all(flags)
          coco_warn(opts, 1, opts.cont.LogLevel, ...
            '%s: %s\n%s: %s\n', mfilename, ...
            'no direction admissible', ...
            'active boundary or terminal constraints were', S.msg);
        end
      else
        [opts atlas cseg] = init_admissible@AtlasBase(atlas, opts, cseg, S);
      end
    end
    
    function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
      chart         = cseg.curr_chart;
      chart.pt_type = 'EP';
      chart.ep_flag = 1;
      if chart.pt<atlas.cont.PtMX
        atlas = atlas.merge_into_atlas({chart});
        flush = false;
      else
        flush = true;
      end
    end
    
    function [prob atlas cseg predict] = refine_step(atlas, prob, cseg)
      [prob atlas cseg predict] = refine_step@AtlasBase(atlas, prob, cseg);
      [atlas chart] = atlas.set_fields(cseg.curr_chart);
      cseg.ptlist   = [ cseg.ptlist atlas.ghost_chart(cseg, chart) ];
    end
    
    function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
      [atlas chart] = atlas.set_fields(cseg.curr_chart);
      [prob cseg]   = cseg.add_chart(prob, chart);
      flush         = true;
      
      if ~atlas.isneighbour(cseg.ptlist{1}, cseg.ptlist{end})
        chart            = atlas.ghost_chart(cseg, cseg.ptlist{end});
        chart.pt_type    = 'GAP';
        chart.ep_flag    = 2;
        cseg.ptlist{end} = chart;
        %cseg.Status      = cseg.CurveSegmentCorrupted;
      end
    end
    
    function [prob atlas cseg correct] = predict(atlas, prob, cseg)
      chart  = atlas.charts{atlas.boundary(1)};
      x0     = chart.x;
      TS     = chart.TS;
      s      = chart.s(:,chart.bv(1));
      h      = atlas.cont.Rmarg*chart.R;
      prcond = struct('x', x0, 'TS', TS, 's', s, 'h', h);
      th     = atlas.cont.theta;
      if th>=0.5 && th<=1
        x1           = x0 + (th*h)*(TS*s);
        [prob cseg]  = CurveSegment.create(prob, chart, prcond, x1);
        [prob ch2]   = cseg.update_TS(prob, cseg.curr_chart);
        s            = ch2.TS'*(x1-x0);
        h            = norm(s);
        s            = s/h;
        x1           = x0 + (h/th)*(ch2.TS*s);
        prcond       = struct('x', x0, 'TS', ch2.TS, 's', s, 'h', h/th);
      else
        x1           = x0 + h*(TS*s);
      end
      [prob chart]   = cseg.update_t(prob, cseg.ptlist{1});
      [prob chart]   = cseg.update_p(prob, chart);
      cseg.ptlist{1} = chart;
      correct        = true;
    end
    
    function [prob atlas cseg accept] = flush(atlas, prob, cseg)
      [atlas last_chart] = atlas.merge_into_atlas(cseg.ptlist(2:end));
      if last_chart && cseg.Status==cseg.CurveSegmentOK
        cseg.ptlist{end}.pt_type = 'EP';
        cseg.ptlist{end}.ep_flag = 1;
      end
      [prob atlas cseg] = atlas.flush@AtlasBase(prob, cseg);
      accept            = isempty(atlas.boundary);
    end
    
  end
  
  methods( Access = private ) % private utility functions
    
    function chart = ghost_chart(atlas, cseg, chart)
      src_chart     = cseg.src_chart;
      s             = src_chart.s(:,src_chart.bv(1));
      h             = atlas.cont.Rmarg*chart.R;
      chart.TS      = src_chart.TS;
      chart.x       = src_chart.x + h*(src_chart.TS*s);
      chart.ep_flag = 2;
    end
    
    function [atlas chart] = set_fields(atlas, chart)
      chart.pt      = atlas.next_pt;
      atlas.next_pt = atlas.next_pt + 1;
      chart.id      = 0;
      chart.R       = atlas.cont.h;
      chart.s       = atlas.cont.s0;
      chart.v       = atlas.cont.v0;
      chart.bv      = atlas.cont.bv0;
      chart.nb      = atlas.cont.nb0;
    end
    
    function flag = isneighbour(atlas, chart1, chart2)
      al  = atlas.cont.almax;
      ta  = tan(al);
      R   = atlas.cont.h;
      x1  = chart1.x;
      x2  = chart2.x;
      dx  = x2-x1;
      x1s = chart1.TS*(chart1.TS'*dx);
      x2s = chart2.TS*(chart2.TS'*dx);
      dst = [ norm(x1s) norm(x2s) norm(dx-x1s) norm(dx-x2s) ...
        subspace( chart1.TS , chart2.TS ) ];
      dstmx = [ R R ta*norm(x1s) ta*norm(x2s) al];
      flag  = all(dst<dstmx);
    end
    
    function flag = isclose(atlas, chart1, chart2)
      al   = atlas.cont.almax;
      R    = atlas.cont.h;
      ta   = tan(al);
      t2a  = tan(2*al);
      x1   = chart1.x;
      x2   = chart2.x;
      dx   = x2-x1;
      phi1 = chart1.TS'*dx;
      phi2 = chart2.TS'*dx;
      x1s  = chart1.TS*(phi1);
      x2s  = chart2.TS*(phi2);
      dst  = [ norm(x1s) norm(x2s) norm(dx-x1s) norm(dx-x2s) ...
        subspace( chart1.TS , chart2.TS ) ];
      n1mx  = ta*min(R,norm(x1s)) + t2a*max(0,norm(x1s)-R);
      n2mx  = ta*min(R,norm(x2s)) + t2a*max(0,norm(x2s)-R);
      dstmx = [ 2*R 2*R n1mx n2mx 2*al];
      if all(dst<dstmx);
        test1 = chart1.v.*(chart1.s'*phi1) - norm(phi1)^2/2;
        test2 = chart2.v.*(chart2.s'*phi2) + norm(phi2)^2/2;
        flag  = any(test1>0) && any(test2<0);
      else
        flag = false;
      end
    end
    
    function [atlas last_chart] = merge_into_atlas(atlas, ptlist)
      last_chart = false;
      bd_charts  = atlas.charts(atlas.boundary);
      for i=1:numel(ptlist)
        chart         = ptlist{i};
        chart.id      = atlas.next_id;
        atlas.next_id = atlas.next_id+1;
        idx           = cellfun(@(x) atlas.isclose(chart, x), bd_charts);
        close_charts  = atlas.boundary(idx);
        checked       = [ 0 chart.id ];
        for k = close_charts
          [atlas chart checked] = atlas.merge_recursive(chart, k, checked);
        end
        atlas.charts   = [ atlas.charts   { chart } ];
        atlas.boundary = [ atlas.boundary  chart.id ];
        bd_charts      = atlas.charts(atlas.boundary);
        idx            = cellfun(@(x) ~isempty(x.bv), bd_charts);
        atlas.boundary = atlas.boundary(idx);
        bd_charts      = bd_charts(idx);
        last_chart     = (chart.pt >= atlas.cont.PtMX);
      end
      if last_chart
        atlas.boundary = [];
      end
      last_chart = isempty(atlas.boundary);
    end
    
    function [atlas chart1 checked] = merge_recursive(atlas, chart1, k, checked)
      if ~any(k==checked)
        checked(end+1) = k;
        chartk = atlas.charts{k};
        if atlas.isclose(chart1, chartk)
          [atlas chart1 chartk] = atlas.voronoi(chart1, chartk, k);
          for k = setdiff(chartk.nb, checked)
            [atlas chart1 checked] = atlas.merge_recursive(chart1, k, checked);
          end
        end
      end
    end
    
    function [atlas chart1 chartk] = voronoi(atlas, chart1, chartk, k)
      dx     = chartk.x - chart1.x;
      phi1   = chart1.TS'*dx;
      phik   = chartk.TS'*(-dx);
      test1  = chart1.v.*(chart1.s'*phi1) - norm(phi1)^2/2;
      testk  = chartk.v.*(chartk.s'*phik) - norm(phik)^2/2;
      flag1  = (test1>0);
      flagk  = (testk>0);
      chart1 = atlas.subtract_half_space(chart1, test1, phi1, flag1, k);
      chartk = atlas.subtract_half_space(chartk, testk, phik, flagk, chart1.id);
      atlas.charts{k} = chartk;
    end
    
    function chart = subtract_half_space(atlas, chart, test, phi, flag, NB)
      k        = find(flag & ~circshift(flag, -1), 1);
      flag     = circshift(flag, -k(1));
      test     = circshift(test, -k(1));
      chart.s  = circshift(chart.s, [0 -k(1)]);
      chart.v  = circshift(chart.v, -k(1));
      chart.nb = circshift(chart.nb, [0 -k(1)]);
      j        = find(~flag & circshift(flag,  -1), 1);
      vx1      = chart.v(j)*chart.s(:,j);
      vx2      = chart.v(j+1)*chart.s(:,j+1);
      nvx1     = vx1 - test(j)/((vx2-vx1)'*phi) * (vx2 - vx1);
      vx1      = chart.v(end)*chart.s(:,end);
      vx2      = chart.v(1)*chart.s(:,1);
      nvx2     = vx1 - test(end)/((vx2 - vx1)'*phi) * (vx2 - vx1);
      chart.s  = [ chart.s(:,1:j) nvx1/norm(nvx1) nvx2/norm(nvx2)];
      chart.v  = [ chart.v(1:j) ; norm(nvx1) ; norm(nvx2)];
      chart.nb = [ chart.nb(1:j+1) NB];
      chart.bv = find( (chart.pt==0 || ~chart.ep_flag) ...
        & (chart.v > atlas.cont.Rmarg*chart.R));
    end
    
  end
  
end