classdef atlas_1d_min < AtlasBase
  %ATLAS_1D_MIN   1d atlas algorithm
  % 
  % This subclass to AtlasBase implements a 1-dimensional expanding
  % boundary atlas algorithm with constant step size and theta-based
  % projection condition and predictor. The algorithm may start from a
  % computational boundary and performs continuation in each direction
  % along the solution manifold until the boundary of the computational
  % domain is reached or the maximum number of steps is reached.
  
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
      % Create an instance of the class object using the class
      % constructor, store class settings in the cont variable, and
      % allocate an empty projection condition.
      
      atlas = atlas_1d_min(prob, cont, dim);
      prob  = CurveSegment.add_prcond(prob, dim);
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
