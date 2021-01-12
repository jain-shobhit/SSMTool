function coco_print_slots(fhan, prob)
%COCO_PRINT_SLOTS([fhan] prob)

if nargin<2
  prob = fhan;
  fhan = 1;
end

if isfield(prob, 'slots')
  t        = coco_opts_tree();
  funcs    = prob.slots.funcs;
  signames = setdiff(fieldnames(prob.slots), 'funcs');
  for i=1:numel(signames)
    signame = signames{i};
    slots   = prob.slots.(signame);
    for idx=slots
      slotid = sprintf('%s@%s', signame, func2str(funcs(idx).F));
      path   = funcs(idx).identifyer;
      if t.prop_exist(path, 'signal', '-no-inherit')
        props = t.prop_get(path, '');
        pr.signals = { props.signal slotid };
        t = t.prop_set(path, '', pr);
      elseif t.prop_exist(path, 'signals', '-no-inherit')
        siglist = t.prop_get(path, 'signals');
        siglist = [ siglist slotid ]; %#ok<AGROW>
        t = t.prop_set(path, 'signals', siglist);
      else
        t = t.prop_set(path, 'signal', slotid);
      end
    end
  end
  t.print_tree(fhan, 'slots', true);
end

end
