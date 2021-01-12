classdef atlas_1d_min < AtlasBase
  %ATLAS_1D_MIN   1d atlas algorithm
  % 
  % This toolbox implements a 1-dimensional expanding boundary atlas
  % algorithm with constant step size and theta-based projection condition
  % and predictor. The algorithm may start from a computational boundary
  % and performs continuation in each direction along the solution manifold
  % until the boundary of the computational domain is reached or the
  % maximum number of steps is reached. The algorithm detects fold- and
  % branch points and uses chart data in order to pass information about
  % the global Jacobian from the linear solver to the atlas algorithm.
  
  % Copyright (C) Frank Schilder, Harry Dankowicz
  % $Id: atlas_1d_min.m 2839 2015-03-05 17:09:01Z fschild $
  
  properties (Access=private)
    boundary = {};       % Cell array to hold atlas boundary charts
    cont     = struct(); % Struct used to hold class settings
  end
  
  methods (Access=private)
    function atlas = atlas_1d_min(prob, cont, dim)
      %ATLAS_1D_MIN   Class constructor.
      %
      % Verify that the continuation problem is initializable with
      % a 1-dimensional atlas algorithm, call the superclass
      % constructor to instantiate an object, and append the cont
      % property with default class settings.
      
      assert(dim==1, '%s: wrong manifold dimension', mfilename);
      atlas      = atlas@AtlasBase(prob);
      atlas.cont = atlas.get_settings(cont);
    end
  end
  
  methods (Static)
    function [prob cont atlas] = create(prob, cont, dim)
      %CREATE   Static construction method.
      %
      % Create an instance of the class object using the class constructor,
      % store class settings in the cont variable, and allocate an empty
      % projection condition. If fold detection is desired, call the
      % add_test_FP constructor after all embeddable monitor functions have
      % been added to the continuation problem structure.
      
      atlas = atlas_1d_min(prob, cont, dim);
      prob  = CurveSegment.add_prcond(prob, dim);
      if atlas.cont.FP
        prob = coco_add_func_after(prob, 'mfunc', ...
          @atlas_1d_min.add_test_FP);
      end
      if atlas.cont.BP
        prob = coco_set(prob, 'lsol', 'det', true); % Store determinant information in chart data
        fid  = coco_get_id('atlas', 'test', 'BP');
        prob = coco_add_func(prob, fid, @atlas.test_BP, [], ...
          'singular', fid, 'uidx', 'all', 'passChart', ... % Depends on all continuation variables and active continuation parameters
          'returnsProb', 'fdim', 1); % Pass chart (to access chart data) and return continuation problem structure
        prob = coco_add_event(prob, 'BP', fid, 0); % BP - Event type
      end
    end
    
    function prob = add_test_FP(prob)
      %ADD_TEST_FP   Append monitor function for fold points.
      %
      % Add monitor function only if there are active continuation
      % parameters. Defer initial call to monitor function by declaring
      % range dimension.
      
      p_idx = coco_get_func_data(prob, 'efunc', 'pidx');
      if numel(p_idx)>=1
        fid  = coco_get_id('atlas', 'test', 'FP');
        prob = coco_add_func(prob, fid, @atlas_1d_min.test_FP, [], ... % Depends on all continuation variables and active continuation parameters
          'singular', fid, 'uidx', 'all', 'passTangent', 'fdim', 1);   % Pass tangent
        prob = coco_add_event(prob, 'FP', fid, 0); % FP - Event type
      end
    end
    
    function [data y] = test_FP(prob, data, u, t)
      %TEST_FP   Monitor function for fold points.
      %
      % Monitor component of tangent vector along primary active
      % continuation parameter.
      
      p_idx = coco_get_func_data(prob, 'efunc', 'pidx');
      y = t(p_idx(1));
    end
    
    function [prob data chart y] = test_BP(prob, data, chart, u)
      %TEST_BP   Monitor function for branch points.
      %
      % Check for singular Jacobian of restricted continuation problem.
      
      cdata = coco_get_chart_data(chart, 'lsol'); % Extract chart data
      if ~isfield(cdata, 'det')
        [prob chart] = prob.cseg.update_det(prob, chart); % Compute determinant, if missing
        cdata = coco_get_chart_data(chart, 'lsol');
      end
      y = cdata.det;
    end
    
  end
  
  methods (Static, Access=private)
    cont = get_settings(cont)
  end
  
  methods (Access=public) % interface methods
    [prob atlas cseg correct] = init_prcond    (atlas, prob, chart)   % Initialize projection condition
    [prob atlas cseg]         = init_chart     (atlas, prob, cseg)    % Initialize chart
    [prob atlas cseg]         = init_admissible(atlas, prob, cseg, S) % Remove inadmissible directions
    [prob atlas cseg flush]   = init_atlas     (atlas, prob, cseg)    % Initialize atlas
    [prob atlas cseg]         = flush          (atlas, prob, cseg)    % Flush point list
    [prob atlas cseg correct] = predict        (atlas, prob, cseg)    % Compute predictor
    [prob atlas cseg flush]   = add_chart      (atlas, prob, cseg)    % Add chart to point list
  end
  
  methods (Access=private)
    flag = isneighbor(atlas, chart1, chart2)
    [atlas cseg] = merge(atlas, cseg)
  end
  
end
