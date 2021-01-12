function coco_print_funcs(fhan, prob)
%COCO_PRINT_FUNCS([fhan] prob)

if nargin<2
  prob = fhan;
  fhan = 1;
end

if isfield(prob, 'efunc') && isfield(prob.efunc, 'funcs')
  t     = coco_opts_tree();
  funcs = prob.efunc.funcs;
  for i=1:numel(funcs)
    func = funcs(i);
    t = t.prop_set(func.identifyer, 'type', func.type);
    if func.pnum>0
      t = t.prop_set(func.identifyer, 'pars', func.pnames);
    end
  end
  t.print_tree(fhan, 'funcs', true);
end

end
