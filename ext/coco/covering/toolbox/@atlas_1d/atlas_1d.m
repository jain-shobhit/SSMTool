classdef atlas_1d < AtlasBase
  
  properties (Access=private)
    PtMX = [], chart_list = {}, func_list = {}, cont = struct();
    first = true;
  end
  
  methods (Access=private) % constructor
    
    function atlas = atlas_1d(prob, cont, dim)
      assert(dim==1, '%s: wrong manifold dimension dim=%d, expected dim=1', ...
        mfilename, dim);
      atlas      = atlas@AtlasBase(prob);
      atlas.cont = atlas.get_settings(cont);
    end
    
  end
  
  methods (Static=true, Access = public)
    [cont, spec] = get_settings(cont)
  end
  
  methods (Static=true) % static construction method
    
    function [prob, cont, atlas] = create(prob, cont, dim)
      atlas      = atlas_1d(prob, cont, dim);
      CurveSegment.class_props('defaults');
      atlas.cont = CurveSegment.class_props(atlas.cont);
      prob       = CurveSegment.add_prcond(prob, dim);
      cont       = atlas.cont;
      
      prob       = coco_add_signal(prob, 'remesh', mfilename);
      
      if atlas.cont.FP
        prob = coco_add_func_after(prob, 'mfunc', @atlas_1d.add_test_FP);
      end
      
      if atlas.cont.BP
        if coco_exist('det', 'class_prop', prob, 'lsol')
          det_flag = coco_get(prob, 'lsol', 'det');
        else
          prob = coco_set(prob, 'lsol', 'det', true);
          det_flag = true;
        end
        if det_flag
          fid  = coco_get_id('atlas', 'test', 'BP');
          prob = coco_add_func(prob, fid, @atlas_1d.test_BP, [], ...
            'singular', fid, 'PassChart', 'fdim', 1);
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
        data.pidx = 1;
        if ~isempty(prob.atlas.cont.fpar)
          data.pidx = find(prob.efunc.acp_idx == ...
            coco_par2idx(prob, prob.atlas.cont.fpar));
        end
        prob = coco_add_func(prob, fid, @atlas_1d.test_FP, ...
          data, 'singular', fid, 'uidx', 'all', 'PassTangent', 'fdim', 1);
        prob = coco_add_event(prob, 'FP', 'SP', fid, 0);
      end
    end
    
    function [data, f] = test_FP(prob, data, x, t) %#ok<INUSL>
      f = t(prob.efunc.p_idx(data.pidx));
    end
    
    function [data, chart, f] = test_BP(prob, data, chart, x) %#ok<INUSD>
      cdata = coco_get_chart_data(chart, 'lsol');
      if ~isfield(cdata, 'det')
        [prob, chart] = prob.cseg.update_det(prob, chart); %#ok<ASGLU>
        cdata = coco_get_chart_data(chart, 'lsol');
      end
      f = cdata.det;
    end
    
  end
  
  methods % interface methods
    [prob, atlas, cseg, correct] = init_prcond    (atlas, prob, chart)
    [prob, atlas, cseg]          = init_chart     (atlas, prob, cseg)
    [prob, atlas, cseg]          = init_admissible(atlas, prob, cseg, S)
    [prob, atlas, cseg, flush]   = init_atlas     (atlas, prob, cseg)
    [prob, atlas, cseg, predict] = refine_step    (atlas, prob, cseg)
    [prob, atlas, cseg, flush]   = add_chart      (atlas, prob, cseg)
    [prob, atlas, cseg, correct] = predict        (atlas, prob, cseg)
    [prob, atlas, cseg]          = flush          (atlas, prob, cseg)
  end
  
end
