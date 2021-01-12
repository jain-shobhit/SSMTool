classdef atlas_1ds < AtlasBase

  properties (Access=private)
    base_chart = struct([]);
    cont       = struct([]);
  end
  
  methods (Access=private) % constructor
    
    function atlas = atlas_1ds(opts, cont, dim)
      assert(dim==1, '%s: wrong manifold dimension dim=%d, expected dim=1', ...
        mfilename, dim);
      atlas      = atlas@AtlasBase(opts);
      atlas.cont = atlas.get_settings(cont);
    end
    
  end
  
  methods (Static=true, Access = private)
    
    function cont = get_settings(cont)
      defaults.PtMX      = 100    ; % number of continuation steps
      defaults.h0        = 0.1    ; % continuation step size, h takes precedence
      defaults.Valpha    = 80     ; % 'vertical' initial tangent
      defaults.FP        = true   ; % detect fold points
      defaults.BP        = true   ; % detect branch points
      defaults.corrector = 'nwtns'; % corrector toolbox
      defaults.interp    = 'cubic'; % use curve segment with cubic interpolation
      cont               = coco_merge(defaults, cont);
      if isfield(cont, 'ItMX')
        cont.PtMX = cont.ItMX;
      end
      assert(any(numel(cont.PtMX)==[1 2]));
      if numel(cont.PtMX) == 2
        if cont.PtMX(1)==0
          cont.PtMX = abs(cont.PtMX(2));
        elseif cont.PtMX(2)==0
          cont.PtMX = -abs(cont.PtMX(1));
        else
          error('%s: only unidirectional continuation is supported', ...
            mfilename);
        end
      end
      if ~isfield(cont, 'h')
        cont.h = cont.h0;
      end
    end
    
  end
    
  methods (Static=true) % static construction method
    
    function [opts cont atlas] = create(opts, cont, dim)
      atlas      = atlas_1ds(opts, cont, dim);
      CurveSegment.class_props('defaults');
      atlas.cont = CurveSegment.class_props(atlas.cont);
      opts       = CurveSegment.add_prcond(opts, dim);
      cont       = atlas.cont;
      
      if atlas.cont.FP
        opts = coco_add_func_after(opts, 'mfunc', @atlas_1ds.add_test_FP);
      end
      
      if atlas.cont.BP
        if coco_exist('det', 'class_prop', opts, 'lsol')
          det_flag = coco_get(opts, 'lsol', 'det');
        else
          opts = coco_set(opts, 'lsol', 'det', true);
          det_flag = true;
        end
        if det_flag
          fid  = coco_get_id('atlas', 'test', 'BP');
          opts = coco_add_func(opts, fid, @atlas_1ds.test_BP, [], ...
            'singular', fid, 'xidx', 'all', 'PassChart', 'fdim', 1);
          opts = coco_add_event(opts, 'BP', 'SP', fid, 0);
        else
          coco_warn(opts, 1, cont.LogLevel, ...
            '%s: %s,\n  %s', mfilename, ...
            'property ''det'' of class ''lsol'' is set to ''false''', ...
            'cannot add BP test function');
        end
      end
      
    end
    
    function opts = add_test_FP(opts)
      if numel(opts.efunc.p_idx)>=1
        fid  = coco_get_id('atlas', 'test', 'FP');
        opts = coco_add_func(opts, fid, @atlas_1ds.test_FP, ...
          [], 'singular', fid, 'xidx', 'all', 'PassTangent', 'fdim', 1);
        opts = coco_add_event(opts, 'FP', 'SP', fid, 0);
      end
    end
    
    function [data f] = test_FP(opts, data, x, t) %#ok<INUSL>
      f = t(opts.efunc.p_idx(1));
    end
    
    function [data chart f] = test_BP(opts, data, chart, x) %#ok<INUSD>
      cdata = coco_get_chart_data(chart, 'lsol');
      if ~isfield(cdata, 'det')
        [opts chart] = opts.cseg.update_det(opts, chart); %#ok<ASGLU>
        cdata = coco_get_chart_data(chart, 'lsol');
      end
      f = cdata.det;
    end
    
  end
  
  methods % interface methods
    
    function [opts atlas cseg correct] = init_prcond(atlas, opts, chart)
      % Initialize initial chart
      chart.R       = 0;
      chart.pt      = -1;
      chart.pt_type = 'IP';
      chart.ep_flag = 1;
      [opts cseg]   = CurveSegment.create_initial(opts, chart, ...
        atlas.cont.Valpha, atlas.cont.h0, 1);
      correct       = cseg.correct;
    end
    
    function [opts atlas cseg flush] = init_atlas(atlas, opts, cseg)
      % Insert chart in point list and flush.
      chart           = cseg.curr_chart;
      chart.s         = sign(atlas.cont.PtMX);
      atlas.cont.PtMX = abs(atlas.cont.PtMX);
      chart.R         = atlas.cont.h;
      chart.pt        = 0;
      chart.pt_type   = 'EP';
      chart.ep_flag   = 1;
      [opts cseg]     = cseg.add_chart(opts, chart);
      flush           = true;
      
      if chart.s~=0
        chart           = cseg.ptlist{1};
        chart.t         = chart.TS*chart.s;
        [opts chart]    = cseg.update_p(opts, chart);
        cseg.ptlist{1}  = chart;
      end
    end
    
    function [opts atlas cseg flush] = add_chart(atlas, opts, cseg)
      chart    = cseg.curr_chart;
      chart.pt = chart.pt + 1;
      if chart.pt >= atlas.cont.PtMX
        chart.pt_type = 'EP';
        chart.ep_flag = 1;
      end
      [opts cseg] = cseg.add_chart(opts, chart);
      flush       = true;
    end
    
    function [opts atlas cseg correct] = predict(atlas, opts, cseg) %#ok<INUSD>
      chart       = atlas.base_chart;
      prcond      = struct('x', chart.x, 'TS', chart.TS, 's', chart.s, 'h', chart.R);
      x1          = chart.x+chart.R*(chart.TS*chart.s);
      [opts cseg] = CurveSegment.create(opts, chart, prcond, x1);
      correct     = true;
    end
    
    function [opts atlas cseg] = flush(atlas, opts, cseg)
      [opts atlas cseg] = atlas.flush@AtlasBase(opts, cseg);
      if cseg.Status == cseg.CurveSegmentOK
        atlas.base_chart = cseg.ptlist{end};
        if atlas.base_chart.pt>=atlas.cont.PtMX
          cseg.Status = cseg.BoundaryPoint;
        end
      end
    end
    
  end
  
end