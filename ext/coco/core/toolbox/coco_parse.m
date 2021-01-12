function [data, opts] = coco_parse(grammar, dspec, ospec, varargin)

str  = coco_stream(varargin{:});
data = struct();
if ~isempty(dspec)
  data = coco_parse_data(grammar, dspec, str);
end
opts = struct();
if ~isempty(ospec)
  opts = coco_parse_opts(ospec, str);
end

end
