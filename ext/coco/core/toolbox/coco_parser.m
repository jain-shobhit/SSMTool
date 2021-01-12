classdef coco_parser < coco_stream
  % COCO_STREAM convert cell array to stream object.
  
  methods
    
    function p = coco_parser(varargin)
      p = p@coco_stream(varargin{:});
      % fprintf(2, '%s: warning: %s\n%s\n', ...
      %   mfilename, 'class coco_parser renamed to coco_stream', ...
      %   'please use coco_stream instead');
    end
    
  end
  
end
