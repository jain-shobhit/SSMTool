classdef coco_cleanup < handle
  
  properties (Access=private)
    fhans = []
    funcs = {}
  end
  
  methods
    
    function key = coco_cleanup()
    end
    
    function fclose(key, fhan)
      key.fhans = [fhan key.fhans];
    end
    
    function call(key, func, varargin)
      key.funcs = [ key.funcs ; {func varargin} ];
    end
    
    function delete(key)
      try %#ok<TRYNC>
        for h=key.fhans
          fclose(h);
        end
        for i=1:size(key.funcs,1)
          func = key.funcs{i,1};
          args = key.funcs{i,2};
          func(args{:});
        end
      end
    end
    
  end
  
end
