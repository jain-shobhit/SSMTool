classdef atlas_0d_recipes < AtlasBase
  %ATLAS_0D_RECIPES   Recipes for continuation: 0d default atlas covering algorithm
  % 
  % This subclass to AtlasBase implements the basic functionality
  % for locating a unique solution to a regular closed continuation
  % problem. Implementation allows for multiple remesh cycles. Each
  % successfully located solution is stored to disk.
  %
  % Class settings: NAdapt    - number of remesh->correct cycles
  % (default: 0) RMMX      - maximum number of remesh loops
  % (default: 10) corrector - nonlinear solver (default: recipes)
  
  % Copyright (C) Frank Schilder, Harry Dankowicz
  % $Id: atlas_0d_recipes.m 2839 2015-03-05 17:09:01Z fschild $
  
  properties
    cont       = struct() % Struct used to hold class settings.
    base_chart = struct() % Struct used to hold a base chart for the predictor.
  end
  
  methods
    
    function atlas = atlas_0d_recipes(opts, cont, dim)
      %ATLAS_0D_RECIPES   Class constructor.
      %
      % Verify that the continuation problem is initializable with
      % a 0-dimensional atlas algorithm, call the superclass
      % constructor to instantiate an object, and append the cont
      % property with default class settings.
      
      assert(dim==0, '%s: wrong manifold dimension dim=%d, expected dim=0', ...
        mfilename, dim);
      atlas      = atlas@AtlasBase(opts);
      atlas.cont = atlas.get_settings(cont);
    end
    
  end
  
  methods (Static=true, Access = private)
    
    function cont = get_settings(cont)
      %GET_SETTINGS   Defines default class settings.
      %
      % Merge default class settings with cont argument.
      
      defaults.NAdapt    =  0; % perform NAdapt remesh->correct cycles
      defaults.RMMX      = 10; % maximum number of remesh loops
      defaults.corrector = 'recipes'; % corrector toolbox
      cont = coco_merge(defaults, cont);
    end
    
  end
  
  methods (Static=true)
    
    function [opts cont atlas] = create(opts, cont, dim)
      %CREATE   Static construction method.
      %
      % Create an instance of the class object using the class
      % constructor, store class settings in the cont variable, and
      % add a remesh signal.
      
      atlas = atlas_0d_recipes(opts, cont, dim);
      cont  = atlas.cont;
      opts  = coco_add_signal(opts, 'remesh', mfilename);
    end
    
  end
  
  methods % interface methods
    
    function [prob atlas cseg correct] = init_prcond(atlas, prob, chart)
      %INIT_PRCOND   Construct an initial curve segment.
      %
      % Create an initial curve segment class instance and proceed
      % to apply a nonlinear corrector. If successful, proceed to
      % init_chart, otherwise proceed to flush with cseg.Status =
      % ~cseg.CurveSegmentOK.
      
      chart.R       = 0;
      chart.pt      = -1;
      chart.pt_type = 'IP';
      chart.ep_flag = 1;
      cseg          = CurveSegmentBase(prob, chart, true);
      correct       = true;
    end
    
    function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
      %INIT_ATLAS   Construct an initial curve segment point list.
      %
      % Function is called following successful execution of
      % init_chart and init_admissible with cseg.Status =
      % cseg.CurveSegmentOK. An initial curve segment point list is
      % created using successfully located initial chart. Execution
      % proceeds to flush.
      
      chart         = cseg.curr_chart;
      chart.pt      = 0;
      chart.pt_type = 'EP';
      chart.ep_flag = 1;
      [prob chart chart.p] = prob.efunc.monitor_F(prob, chart, chart.x, chart.t);
      cseg.ptlist   = { chart };
      flush         = true;
    end
    
    function [prob atlas cseg] = flush(atlas, prob, cseg)
      %FLUSH   Write point list to screen and disk.
      %
      % Call superclass flush function to write point list
      % (consisting of single chart) to screen and disk. If point
      % number exceeds NAdapt, then terminate. Otherwise proceed to
      % predict.
      
      [prob atlas cseg] = atlas.flush@AtlasBase(prob, cseg, 'all');
      if cseg.Status == cseg.CurveSegmentOK
        atlas.base_chart = cseg.ptlist{end};
        if atlas.base_chart.pt>=atlas.cont.NAdapt
          cseg.Status = cseg.BoundaryPoint;
        end
      end
    end
    
    function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
      %ADD_CHART   Construct a curve segment point list.
      %
      % A curve segment point list is created using successfully
      % located initial chart. Execution proceeds to flush.
      
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
    
    function [opts atlas cseg correct] = predict(atlas, opts, cseg) %#ok<INUSD>
      %PREDICT   Construct predictor.
      %
      % Apply the remesh functions associated with the embedded
      % continuation problem objects to the base chart and
      % construct new curve segment. Proceed to apply a nonlinear
      % corrector. If successful, proceed to add_chart, otherwise
      % proceed to flush with cseg.Status = ~cseg.CurveSegmentOK.
      
      chart = atlas.base_chart;
      
      x0              = chart.x;
      [opts chart x0] = coco_remesh(opts, chart, x0, [], atlas.cont.RMMX);
      chart.x         = x0;
      
      cseg    = CurveSegmentBase(opts, chart, false);
      correct = true;
    end
    
  end
  
end
