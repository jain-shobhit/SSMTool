classdef atlas_kd < AtlasBase
  %ATLAS_KD   k-d atlas algorithm
  %
  % This subclass to AtlasBase implements a k-dimensional expanding
  % boundary atlas algorithm with adaptive step size and theta-based
  % projection condition and predictor that recognizes computational domain
  % boundaries and may start from a boundary. The algorithm is designed to
  % prevent redundant coverage and premature termination.
  
  % This file is part of the atlas_kd toolbox, copyright (C) Michael
  % Henderson, Frank Schilder, Harry Dankowicz, Erika Fotsch, of the package
  % COCO (http://sourceforge.net/projects/cocotools).
  
  properties (Access=private)
    dim      = 0;        % Manifold dimension
    tree     = [];       % The hierarchical bounding box
    boundary = [];       % Array of ids of boundary charts
    charts   = {};       % Cell array of charts in this atlas
    edge     = 0;        % Edge counter
    cont     = struct(); % Struct used to hold class settings
    next_pt  = 0;        % Chart counter
    first    = true;
    funcs    = {};
  end
  
  methods (Access=private)
    function atlas = atlas_kd(prob, cont, dim)
      %ATLAS_KD   Class constructor.
      %
      % Call the superclass constructor to instantiate an object, and
      % append the cont property with default class settings.
      
      atlas      = atlas@AtlasBase(prob);
      atlas.dim  = dim;
      atlas.edge = 2*dim;
      atlas.cont = atlas.get_settings(cont, dim);
    end
  end
  
  methods (Static)
    function [prob, cont, atlas] = create(prob, cont, dim)
      %CREATE   Static construction method.
      %
      % Create an instance of the class object using the class
      % constructor, store class settings in the cont variable, and
      % allocate an empty projection condition. Append a slot function to
      % store atlas to disk.
      
      atlas = atlas_kd(prob, cont, dim);
      prob  = CurveSegment.add_prcond(prob, dim);
      prob  = coco_add_slot(prob, 'atlas', @atlas.save_atlas, ...
        [], 'save_bd');
      prob  = coco_add_signal(prob, 'remesh', mfilename);

    end
    
    function [data, res] = save_atlas(prob, data, varargin)
      %SAVE_ATLAS   Slot function: Save atlas for plotting.
      %
      % Store the charts and boundary properties with the 'atlas'
      % identifier in the bd.mat file.
      
      res.charts   = prob.atlas.charts;
      res.boundary = prob.atlas.boundary;
    end
  end
  
  methods (Static=true, Access = public)
    [cont, spec] = get_settings(cont, dim)
  end
  
  methods (Access=public) % interface methods
    [prob, atlas, cseg, correct] = init_prcond    (atlas, prob, chart)   % Initialize projection condition
    [prob, atlas, cseg]          = init_chart     (atlas, prob, cseg)    % Initialize chart
    [prob, atlas, cseg]          = init_admissible(atlas, prob, cseg, S) % Remove inadmissible directions
    [prob, atlas, cseg, flush]   = init_atlas     (atlas, prob, cseg)    % Initialize atlas
    [prob, atlas, cseg]          = flush          (atlas, prob, cseg)    % Flush point list
    [prob, atlas, cseg, correct] = predict        (atlas, prob, cseg)    % Compute predictor
    [prob, atlas, cseg, flush]   = add_chart      (atlas, prob, cseg)    % Add chart to point list
  end
  
  methods (Access=public)
    atlas = update_atlas(atlas, ID)
    atlas = add_chart2tree(atlas, chart)
    flag  = isneighbor(atlas, chart1, chart2)
    flag  = isclose(atlas, chart1, chart2)
    [prob, atlas, cseg] = merge(atlas, prob, cseg)
  end
    
end
