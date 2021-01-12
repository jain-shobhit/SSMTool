classdef continex_shape
  
  methods (Access=public,Abstract=true)
    f     = F         (shape, u)
    J     = DFDU      (shape, u)
    shape = update    (shape, prcond)
    u     = pull_back (shape, u)
  end
  
  properties (Access=public)
    t = [], u1 = []
    scale = 1, h = 0, MaxAbsDist = inf
  end
  
  methods (Access=public)
    
    function shape = continex_shape(cont) %#ok<INUSD>
    end
    
    function shape = set_scale(shape, s)
      shape.scale = s;
    end
    
  end
  
  methods (Access=public)
    
    function plot_shape(varargin)
    end
    
    function plot_point(varargin)
    end
    
  end
  
  methods (Access=public,Static=true)
    function cont = get_settings(cont)
    end
  end
  
end
