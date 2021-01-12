function varargout = coco_read_tb_info(oid, run, varargin)
%COCO_READ_TB_INFO   Read toolbox info from data file.
%
% VARARGOUT = COCO_READ_TB_INFO(OID, RUN, [LAB], [PROP] ...)
% Read toolbox information from bd-file of run RUN. If LAB is given, read
% toolbox information from solution file of solution LAB of run RUN. If
% names of properties are given with PROP ..., return only values of these
% properties in the order requested. If a property does not exist, an empty
% string is returned.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_read_tb_info.m 3023 2017-08-22 21:28:21Z hdankowicz $

fid = coco_get_id(oid, 'tb_info');
s   = coco_stream(varargin{:});
if isnumeric(s.peek())
  lab     = s.get();
  tb_info = coco_read_solution(fid, run, lab, 'data');
else
  tb_info = coco_bd_read(run, fid);
end

if isempty(s)
  varargout = { tb_info };
else
  nout = 1;
  while ~isempty(s)
    name = s.get();
    if isfield(tb_info, name)
      val = tb_info.(name);
    else
      val = '';
    end
    varargout{nout} = val;
    nout = nout+1;
  end
end

end
