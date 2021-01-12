classdef coco_func_data_ptr < handle
  % Simple pointer class used by coco_func_data.
  
  properties (Access=public)
    pr=struct()
    sh=struct()
  end
  
  methods
    function obj = coco_func_data_ptr(p,s)
      if nargin>0
        obj.pr = p;
        obj.sh = s;
      end
    end
  end
  
end
