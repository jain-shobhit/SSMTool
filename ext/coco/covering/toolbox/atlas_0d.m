classdef atlas_0d < AtlasBase
  
  properties
    cont       = struct()
    base_chart = struct()
  end
  
  methods % constructor
    
    function atlas = atlas_0d(opts, cont, dim)
      assert(dim==0, '%s: wrong manifold dimension dim=%d, expected dim=0', ...
        mfilename, dim);
      atlas      = atlas@AtlasBase(opts);
      atlas.cont = atlas.get_settings(cont);
    end
    
  end
  
  methods (Static=true, Access = private)
    
    function cont = get_settings(cont)
      defaults.NAdapt    =  0; % perform NAdapt remesh->correct cycles
      defaults.RMMX      = 10; % maximum number of remesh loops
      cont = coco_merge(defaults, cont);
    end
    
  end
  
  methods (Static=true) % static construction method
    
    function [opts cont atlas] = create(opts, cont, dim)
      atlas = atlas_0d(opts, cont, dim);
      cont  = atlas.cont;
      opts  = coco_add_signal(opts, 'remesh', mfilename);
    end
    
  end
  
  methods % interface methods
    
    function [prob atlas cseg correct] = init_prcond(atlas, prob, chart)
      chart.R       = 0;
      chart.pt      = -1;
      chart.pt_type = 'IP';
      chart.ep_flag = 1;
      cseg          = CurveSegmentBase(prob, chart, true);
      correct       = true;
      prob = coco_add_signal(prob, 'update', mfilename); % Harry added
    end
    
    function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
      chart         = cseg.curr_chart;
      chart.pt      = 0;
      chart.pt_type = 'EP';
      chart.ep_flag = 1;
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t);
      cseg.ptlist   = { chart };
      flush         = true;
    end
    
    function [prob atlas cseg] = flush(atlas, prob, cseg)
      [prob atlas cseg] = atlas.flush@AtlasBase(prob, cseg, 'all');
      if cseg.Status == cseg.CurveSegmentOK
        atlas.base_chart = cseg.ptlist{end};
        if atlas.base_chart.pt>=atlas.cont.NAdapt
          cseg.Status = cseg.BoundaryPoint;
        end
      end
    end
    
    function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
      chart    = cseg.curr_chart;
      chart.pt = chart.pt + 1;
      if chart.pt >= atlas.cont.NAdapt
        chart.pt_type = 'EP';
        chart.ep_flag = 1;
      end
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t);
      cseg.ptlist = { chart };
      flush       = true;
    end
    
    function [opts atlas cseg correct] = predict(atlas, opts, cseg)
      chart = atlas.base_chart;
      
      % remesh base chart of new curve segment
      x0              = chart.x;
      [opts chart x0] = coco_remesh(opts, chart, x0, [], atlas.cont.RMMX);
      chart.x         = x0;
      
      % construct new curve segment
      cseg    = CurveSegmentBase(opts, chart, false);
      opts    = coco_emit(opts, 'update', cseg); % Harry added
      correct = true;
    end
    
  end
  
end
