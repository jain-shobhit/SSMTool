classdef coco_ptr < coco_func_data
  % COCO_PTR alias class for coco_func_data.
  
  methods
    function obj = coco_ptr(varargin)
      obj = obj@coco_func_data(varargin{:});
    end
  end
end
