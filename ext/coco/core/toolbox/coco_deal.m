function varargout = coco_deal(varargin)
if isa(varargin{1}, 'coco_stream')
  s = coco_stream(varargin{:});
  [varargout{1:nargout}] = s.get;
else
  varargout = varargin(1:nargout);
end
end
