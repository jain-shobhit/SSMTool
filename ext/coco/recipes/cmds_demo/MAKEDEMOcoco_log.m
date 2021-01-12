coco_log(prob, 1, 1, '%s: remeshed, N=%d, H=%d, status=''%s''\n', mfilename, N, H, stat);
[data Jt] = coco_ezDFDX('f(o,d,x)', prob, data, @calcvar_F, u);