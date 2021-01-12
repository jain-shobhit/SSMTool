function prob_new = coco_prob(prob)
%COCO_PROB([prob]) create new or copy of coco problem.

if nargin<1
  prob_new = coco_set();
else
  prob_new = coco_set(prob);
  prob_new = coco_save_funcs(prob_new);
end

end
