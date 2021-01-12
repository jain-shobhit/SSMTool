classdef atlas_1d_recipes < AtlasBase
  %ATLAS_1D_RECIPES   Recipes for continuation: 1d default atlas covering algorithm
  % 
  % This subclass to AtlasBase implements the basic functionality
  % for locating a unique solution to a regular closed continuation
  % problem. Implementation allows for multiple remesh cycles.
  % Each successfully located solution is stored to disk.
  %
  % Class settings:
  % PtMX      - maximum number of continuation steps (default: 50)
  % bi_direct - boolean flag for bidirectional continuation (default: true)
  % interp    - identifier for interpolation scheme in cseg (default: 'cubic') 
  % Valpha    - angle for 'vertical' initial tangent (default: 80)
  % h0        - initial continuation step size (default: 0.1)
  % h_max     - maximum continuation step size (default: 0.5)
  % h_min     - minimum continuation step size (default: 0.01)
  % h_fax max - maximum step size adaptation factor (default: 2.0)
  % h_fac_min - minimum step size adaptation factor (default: 0.5)
  % MaxRes    - maximal residual in prediction step (default: 0.1)
  % al_max    - maximum angle between consecutive tangents (default: 7.0)
  % ga        - adaptation security factor (default: 0.95)
  % FP        - boolean flag for fold point detection (default: false)
  % BP        - boolean flag for branch point detection (default: false)
  % NAdapt    - remesh frequency, 0 == off (default: 0)
  % RMMX      - maximum number of remesh loops (default: 10)
  % corrector - nonlinear solver (default: recipes)
  %
  % Continuation parameters:
  % atlas.test.FP - optional nonembedded continuation parameter associated with detection of fold points.
  % atlas.test.BP - optional nonembedded continuation parameter associated with detection of branch points.

  % Copyright (C) Frank Schilder, Harry Dankowicz
  % $Id: atlas_1d_recipes.m 2839 2015-03-05 17:09:01Z fschild $
  
  properties (Access=private)
    PtMX = [] % Maximum number of continuation steps for each desired direction.
    chart_list = {} % Cell array used to hold a last of active atlas charts.
    func_list = {}
    cont = struct() % Struct used to hold class settings.
  end
  
  methods (Access=private) % constructor
    
    function atlas = atlas_1d_recipes(prob, cont, dim)
      %ATLAS_1D_RECIPES   Class constructor.
      %
      % Verify that the continuation problem is initializable with
      % a 1-dimensional atlas algorithm, call the superclass
      % constructor to instantiate an object, and append the cont
      % property with default class settings.

      assert(dim==1, '%s: wrong manifold dimension dim=%d, expected dim=1', ...
        mfilename, dim);
      atlas      = atlas@AtlasBase(prob);
      atlas.cont = atlas.get_settings(cont);
    end
    
  end
  
  methods (Static=true, Access = private)
    
    function cont = get_settings(cont)
      %GET_SETTINGS   Defines default class settings.
      %
      % Merge default class settings with cont argument.

      defaults.PtMX      = 50     ; % max. number of continuation steps
      defaults.bi_direct = true   ; % go in both directions by default
      defaults.interp    = 'cubic'; % cubic interpolation in cseg
      defaults.Valpha    = 80     ; % 'vertical' initial tangent
      defaults.h0        = 0.1    ; % initial continuation step size
      defaults.h_max     = 0.5    ; % max. continuation step size
      defaults.h_min     = 0.01   ; % min. continuation step size
      defaults.h_fac_max = 2.0    ; % max. step size adaption factor
      defaults.h_fac_min = 0.5    ; % min. step size adaption factor
      defaults.MaxRes    = 0.1    ; % max. residuum for prediction step
      defaults.al_max    = 7.0    ; % max. angle between two consecutive tangents
      defaults.ga        = 0.95   ; % adaption security factor
      defaults.FP        = false  ; % detect fold points
      defaults.BP        = false  ; % detect branch points
      defaults.NAdapt    = 0      ; % adapt every NAdapt steps, 0 == off
      defaults.RMMX      = 10     ; % maximum number of remesh loops
      defaults.corrector = 'recipes'; % corrector toolbox
      
      cont = coco_merge(defaults, cont);
      cont.h0        = min(max(cont.h0, cont.h_min), cont.h_max);
      cont.arc_alpha = cont.al_max * pi / 180;
      if isfield(cont, 'ItMX')
        cont.PtMX = cont.ItMX;
      end
      % Reconstruct a two-element PtMX vector corresponding to
      % maximum number of continuation steps in each direction.
      assert(any(numel(cont.PtMX)==[1 2]));
      if numel(cont.PtMX)==1
        if cont.bi_direct
          cont.PtMX = [ cont.PtMX cont.PtMX ];
        else
          if cont.PtMX>=0
            cont.PtMX = [ 0 cont.PtMX ];
          else
            cont.PtMX = [ cont.PtMX 0 ];
          end
        end
      end
      cont.PtMX = [-1 1] .* abs(cont.PtMX);
    end
    
  end
  
  methods (Static=true)
    
    function [prob cont atlas] = create(prob, cont, dim)
      %CREATE   Static construction method.
      %
      % Create an instance of the class object using the class
      % constructor, store class settings in the cont variable,
      % append an empty projection condition to the restricted
      % continuation problem, and add a remesh signal. Optionally
      % add test functions and event detection for fold and branch
      % points.

      atlas      = atlas_1d_recipes(prob, cont, dim);
      CurveSegment.class_props('defaults');
      atlas.cont = CurveSegment.class_props(atlas.cont);
      prob       = CurveSegment.add_prcond(prob, dim);
      cont       = atlas.cont;
      
      prob       = coco_add_signal(prob, 'remesh', mfilename);
      
      % Use coco_add_func_after to delay construction until after
      % all monitor functions have been appended (at which point
      % prob.efunc.p_idx has been assigned).
      if atlas.cont.FP
        prob = coco_add_func_after(prob, 'mfunc', @atlas_1d_recipes.add_test_FP);
      end
      
      % Branch point detection corresponds to a singular point
      if atlas.cont.BP
        if coco_exist('det', 'class_prop', prob, 'lsol')
          det_flag = coco_get(prob, 'lsol', 'det');
        else
          prob = coco_set(prob, 'lsol', 'det', true);
          det_flag = true;
        end
        if det_flag
          fid  = coco_get_id('atlas', 'test', 'BP');
          prob = coco_add_func(prob, fid, @atlas_1d_recipes.test_BP, ...
            atlas.cont, 'singular', fid, 'xidx', 'all', 'PassChart', 'fdim', 1);
          prob = coco_add_event(prob, 'BP', 'SP', fid, 0);
        else
          coco_warn(prob, 1, cont.LogLevel, ...
            '%s: %s,\n  %s', mfilename, ...
            'property ''det'' of class ''lsol'' is set to ''false''', ...
            'cannot add BP test function');
        end
      end
    end
    
    function prob = add_test_FP(prob)
      if numel(prob.efunc.p_idx)>=1
        fid  = coco_get_id('atlas', 'test', 'FP');
        prob = coco_add_func(prob, fid, @atlas_1d_recipes.test_FP, ...
          [], 'singular', fid, 'xidx', 'all', 'PassTangent', 'fdim', 1);
        prob = coco_add_event(prob, 'FP', 'SP', fid, 0);
      end
    end
    
    function [data f] = test_FP(prob, data, x, t) %#ok<INUSL>
      %TEST_FP   Monitor component of tangent vector.
      % 
      % Returns component of tangent vector along first active
      % continuation parameter.
      
      f = t(prob.efunc.p_idx(1));
    end
    
    function [data chart f] = test_BP(prob, data, chart, x) %#ok<INUSD>
      %TEST_BP   Monitor rescaled determinant of Jacobian.
      %
      % Rescaled determinant is stored in chart data by the
      % lsol_recipes linear solver.
      
      cdata = coco_get_chart_data(chart, 'lsol');
      if isfield(cdata, 'det')
        f = cdata.det;
      else %bug: function does not return prob!
        [prob chart] = prob.cseg.update_det(prob, chart); %#ok<ASGLU>
        cdata = coco_get_chart_data(chart, 'lsol');
        f = cdata.det;
      end
    end
    
  end
  
  methods % interface methods
    
    function [prob atlas cseg correct] = init_prcond(atlas, prob, chart)
      %INIT_PRCOND   Construct an initial curve segment.
      %
      % Create an initial curve segment class instance, with
      % corresponding projection condition, and proceed to apply a
      % nonlinear corrector if all components of chart.t equal 0.
      % If successful, proceed to init_chart, otherwise proceed to
      % flush with cseg.Status = ~cseg.CurveSegmentOK.

      chart.R       = 0;
      chart.pt      = -1;
      chart.pt_type = 'IP';
      chart.ep_flag = 1;
      [prob cseg]   = CurveSegment.create_initial(prob, chart, ...
        atlas.cont.Valpha, atlas.cont.h0, 0);
      correct       = cseg.correct;
    end
    
    function [prob atlas cseg] = init_chart(atlas, prob, cseg)
      %INIT_CHART   Initialize initial chart.
      %
      % Allow for bidirectional continuation. For each available
      % direction, store corresponding insertion method and maximum
      % number of continuation steps. Append initial chart to curve
      % segment point list and log angle of tangent vector with
      % coordinate axes of each active continuation parameter.
      
      chart      = cseg.curr_chart;
      chart.pt   = 0;
      chart.s    = [-1 1];
      chart.im   = {'prepend' 'append'};
      atlas.PtMX = abs(atlas.cont.PtMX);
      if sum(atlas.PtMX)==0
        chart.s    = 1;
        chart.im   = {'append'};
        atlas.PtMX = 0;
      elseif any(atlas.PtMX==0)
        idx        = atlas.PtMX~=0;
        chart.s    = chart.s(idx);
        chart.im   = chart.im(idx);
        atlas.PtMX = atlas.PtMX(idx);
      end
      [prob cseg]     = cseg.add_chart(prob, chart);
      cseg.curr_chart = cseg.ptlist{1};
      
      CurveSegment.log_angle(prob, 2, atlas.cont.LogLevel, ...
        cseg.curr_chart.t, prob.efunc.p_idx, mfilename, 'init_chart');
    end
    
    function [prob atlas cseg] = init_admissible(atlas, prob, cseg, S)
      %INIT_ADMISSIBLE   Compute subset of admissible directions
      %
      % Restrict attention to initial points inside the
      % computational domain or points on the domain boundary with
      % admissible directions that point into the computational
      % domain and proceed to init_atlas. If initial point lies
      % outside of computational domain, the superclass
      % init_admissible function sets cseg.status to
      % CurveSegmentCorrupted and proceeds to flush.
      
      if S.ep_flag == 1
        flags = (S.dir_flags==atlas.IsAdmissible);
        cseg.curr_chart.s  = cseg.curr_chart.s(flags);
        cseg.curr_chart.im = cseg.curr_chart.im(flags);
        atlas.PtMX         = atlas.PtMX(flags);
        if ~any(flags)
          coco_warn(prob, 1, prob.cont.LogLevel, ...
            '%s: %s\n%s: %s\n', mfilename, ...
            'no direction admissible', ...
            'active boundary or terminal constraints were', S.msg);
        end
      else
        [prob atlas cseg] = init_admissible@AtlasBase(atlas, prob, cseg, S);
      end
    end
    
    function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
      %INIT_ATLAS   Initialize atlas chart list
      %
      % For each allowable direction, assign suitable restricted
      % copy of initial chart to atlas chart list and proceed to flush.
      
      chart          = cseg.curr_chart;
      chart.R        = atlas.cont.h0;
      chart.pt_type  = 'EP';
      chart.ep_flag  = 1;
      cseg.ptlist{1} = chart;
      
      % split chart into chart list
      for i=1:numel(chart.s)
        chart1              = chart;
        chart1.s            = chart.s(i);
        chart1.im           = chart.im{i};
        chart1.t            = chart.s(i)*chart.t;
        [prob chart1]       = cseg.update_p(prob, chart1);
        chart1.ignore_at    = (i>1);
        atlas.chart_list{i} = chart1;
        if i==1
          cseg.ptlist{1} = atlas.chart_list{1};
          prob = AtlasBase.bddat_set(prob, 'ins_mode', atlas.chart_list{1}.im);
        else
          [prob atlas.func_list{i}] = coco_save_funcs(prob);
        end
      end
      
      % Flush curve segment.
      flush = true;
    end
    
    function [opts atlas cseg predict] = refine_step(atlas, opts, cseg)
      % correction failed, refine step size
      chart = cseg.ptlist{1};
      if chart.R <= atlas.cont.h_min
        % no convergence with minimum step size -> stop continuation in
        % current direction
        predict = false;
      else
        % retry with reduced step size
        chart.R = max(atlas.cont.h_min, atlas.cont.h_fac_min*chart.R);
        atlas.chart_list{1} = chart;
        predict = true;
        coco_warn(opts, 3, atlas.cont.LogLevel, ...
          'no convergence, retry with step size h = % .2e\n', chart.R);
      end
    end
    
    function [opts atlas cseg flush] = add_chart(atlas, opts, cseg)
      chart    = cseg.curr_chart;
      chart.pt = chart.pt + 1;
      if chart.pt >= atlas.PtMX(1)
        chart.pt_type = 'EP';
        chart.ep_flag = 1;
      end
      [opts cseg]  = cseg.add_chart(opts, chart);
      flush        = true;
      
      cont = atlas.cont; %#ok<PROP>
      
      % check angle condition and refine step size if necessary
      % This test is a modification of `arc_beta > opts.cont.arc_alpha' and
      % allows for accepting charts with angles somewhat larger than arc_alpha.
      % This modification significantly reduces the amount of unnecessarily
      % rejected steps. The 'abs' is necessary to avoid complex numbers
      % that occur sometimes due to roundoff errors.
      chart1   = cseg.ptlist{1};
      chart2   = cseg.ptlist{2};
      t1       = chart1.t;
      R1       = chart1.R;
      arc_beta = abs(acos(t1' * chart2.t));
      coco_log(opts, 2, cont.LogLevel, '%s: angle(t1,t0) = %.2e[DEG]\n', ...
        mfilename, 180 * arc_beta / pi); %#ok<PROP>
      if (chart2.pt>0) && (arc_beta > cont.h_fac_max * cont.arc_alpha) %#ok<PROP>
        
        % reduce step size if possible and repeat continuation step
        if R1 > cont.h_min %#ok<PROP>
          coco_warn(opts, 3, cont.LogLevel, ...
            'atlas: beta [%.4e] > h_fac_max * al_max, refining step size\n', ...
            180 * arc_beta / pi); %#ok<PROP>
          
          chart1.R            = max(cont.h_min, cont.h_fac_min*R1); %#ok<PROP>
          atlas.chart_list{1} = chart1;
          flush               = false;
          return
          
        else % R1 <= cont.h_min
          coco_warn(opts, 3, cont.LogLevel, ...
            'atlas: minimum stepsize reached, but beta [%.4e] > h_fac_max * al_max\n', ...
            180 * arc_beta / pi); %#ok<PROP>
        end
        
      end
      
      % copmpute new radius
      if chart2.pt>0
        % This if-statement takes care of the case arc_beta==0, which produced an
        % error in some versions of Matlab.
        if cont.h_fac_max^2 * arc_beta < cont.arc_alpha %#ok<PROP>
          h_fac = cont.h_fac_max; %#ok<PROP>
        else
          h_fac = cont.arc_alpha / (sqrt(cont.h_fac_max) * arc_beta); %#ok<PROP>
          h_fac = max(cont.h_fac_min, min(cont.h_fac_max, h_fac)); %#ok<PROP>
        end
        chart2.R = cont.ga * h_fac * R1; %#ok<PROP>
        chart2.R = max( min(chart2.R, cont.h_max), cont.h_min); %#ok<PROP>
      end
      
      % check residuum at new vertex
      opts = coco_emit(opts, 'set_mode', 'check_res');
      while chart2.R>cont.h_min %#ok<PROP>
        x = chart2.x + chart2.TS*(chart2.R*chart2.s);
        [opts chart2 f] = opts.efunc.F(opts, chart2, x);
        if norm(f) > cont.MaxRes %#ok<PROP>
          coco_warn(opts, 3, cont.LogLevel, ...
            'atlas: norm(f)=%.4e, refining step size\n', norm(f)); %#ok<PROP>
          h_fac   = max(cont.h_fac_min, cont.MaxRes/norm(f)); %#ok<PROP>
          chart2.R = max(cont.h_min, h_fac*chart2.R); %#ok<PROP>
        else
          break;
        end
      end
      
      % copy chart back as end-point of ptlist
      cseg.ptlist{2} = chart2;
    end
    
    function [opts atlas cseg correct] = predict(atlas, opts, cseg)
      chart = atlas.chart_list{1};
      
      % remesh base chart of new curve segment
      log_angle = [];
      if isempty(cseg)
        nad = atlas.cont.NAdapt;
        if nad > 0 && mod(chart.pt,nad) == 0
          x0                 = chart.x;
          V0                 = [ chart.t chart.TS ];
          [opts chart x0 V0] = coco_remesh(opts, chart, x0, V0, atlas.cont.RMMX);
          nv                 = repmat(sqrt(sum(V0.^2,1)), size(V0,1), 1);
          chart.x            = x0;
          chart.t            = V0(:,1)./nv(:,1);
          chart.TS           = V0(:,2:end)./nv(:,2:end);
          log_angle = coco_log(opts, 3, atlas.cont.LogLevel);
        end
      end
      
      % use pseudo arc-length projection condition
      prcond.x  = chart.x;
      prcond.TS = chart.TS;
      prcond.s  = chart.s;
      prcond.h  = chart.R;
      x1        = chart.x + chart.TS*(chart.R*chart.s);
      % construct new curve segment
      [opts cseg] = CurveSegment.create(opts, chart, prcond, x1);
      correct     = true;
      
      % compute and print angle between remeshed and actual tangent
      if ~isempty(log_angle)
        t0           = chart.t;
        [opts chart] = cseg.update_TS(opts, chart);
        [opts chart] = cseg.update_t(opts, chart);
        t1           = chart.t;
        fprintf(log_angle, '%s: angle(t_remesh,t) = %.2e[DEG]\n', ...
          mfilename, 180*subspace(t0,t1)/pi);
      end
    end
    
    function [opts atlas cseg] = flush(atlas, opts, cseg)
      %FLUSH   Flush chart to disk and screen output.
      
      % flush point list
      [opts atlas cseg] = atlas.flush@AtlasBase(opts, cseg);
      
      if cseg.Status == cseg.CurveSegmentOK
        chart = cseg.ptlist{end};
        
        % flush last point into chart_list
        atlas.chart_list{1} = chart;
        
        % end of branch if It>=PtMX.
        if isempty(atlas.PtMX) || (chart.pt>=atlas.PtMX(1))
          cseg.Status = cseg.BoundaryPoint;
        end
      end
      
      % Handle end of branch.
      if cseg.Status~=cseg.CurveSegmentOK && numel(atlas.chart_list)>1
        atlas.PtMX       = atlas.PtMX(2:end);
        atlas.chart_list = atlas.chart_list(2:end);
        chart            = atlas.chart_list{1};
        atlas.func_list  = atlas.func_list(2:end);
        opts = coco_restore_funcs(opts, atlas.func_list{1});
        opts = AtlasBase.bddat_set(opts, 'ins_mode', chart.im);
        atlas.PrintHeadLine = true;
        cseg.Status=cseg.CurveSegmentOK;
      end
      
    end
    
  end
  
end
