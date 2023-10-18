function prob = cocoSet(opts, prob)
%COCOSET This function present settings for coco

default = cocoOptions();
% settings for continuation
if opts.NPR~=default.NPR
    prob = coco_set(prob, 'cont', 'NPR', opts.NPR);
end
if opts.NSV~=default.NSV
    prob = coco_set(prob, 'cont', 'NSV', opts.NSV);
end
if opts.NAdapt~=default.NAdapt
    prob = coco_set(prob, 'cont', 'NAdapt', opts.NAdapt);
end
if opts.h0~=default.h0
    prob = coco_set(prob, 'cont', 'h0', opts.h0);
end
if opts.h_max~=default.h_max
    prob = coco_set(prob, 'cont', 'h_max', opts.h_max);
end
if opts.h_min~=default.h_min
    prob = coco_set(prob, 'cont', 'h_min', opts.h_min);
end
if opts.h_fac_min~=default.h_fac_min
    prob = coco_set(prob, 'cont', 'h_fac_min', opts.h_fac_min);
end
if opts.MaxRes~=default.MaxRes
    prob = coco_set(prob, 'cont', 'MaxRes', opts.MaxRes);
end
if opts.bi_direct~=default.bi_direct
    prob = coco_set(prob, 'cont', 'bi_direct', opts.bi_direct);
end
if any(opts.PtMX~=default.PtMX)
    prob = coco_set(prob, 'cont', 'PtMX', opts.PtMX);
end
% settings for correction
if opts.ItMX~=default.ItMX
    prob = coco_set(prob, 'corr', 'ItMX', opts.ItMX);
end
if opts.TOL~=default.TOL
    prob = coco_set(prob, 'corr', 'TOL', opts.TOL);
end
% settings for collocation
if opts.NTST~=default.NTST
    prob = coco_set(prob, 'coll', 'NTST', opts.NTST);
end
if opts.NCOL~=default.NCOL
    prob = coco_set(prob, 'coll', 'NCOL', opts.NCOL);
end
if opts.MXCL~=default.MXCL
    prob = coco_set(prob, 'coll', 'MXCL', opts.MXCL);
end
% settings for kd continuation
if opts.atlasdim>1
    if opts.MaxRes~=default.MaxRes
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'MaxRes', opts.MaxRes);
    end
    if opts.PtMX~=default.PtMX
        prob = coco_set(prob, 'cont', 'atlas', 'kd', 'PtMX', opts.PtMX);
    end
    if opts.R~=default.R
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'R', opts.R);
    end
    if opts.R_max~=default.R_max
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'R_max', opts.R_max);
    end
    if opts.R_min~=default.R_min
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'R_min', opts.R_min);
    end
    if opts.R_fac_max~=default.R_fac_max
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'R_fac_max', opts.R_fac_max);
    end
    if opts.R_fac_min~=default.R_fac_min
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'R_fac_min', opts.R_fac_min);
    end
    if opts.ga~=default.ga
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'ga', opts.ga);
    end
    if opts.almax~=default.almax
        prob = coco_set(prob, 'cont', 'atlas', 'kd',  'almax', opts.almax);
    end
end
end