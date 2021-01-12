function coco_print_sigs(varargin)
%COCO_PRINT_SIGS([fhan] prob ['all'])

fhan = 1;
all = false;
s = coco_stream(varargin{:});
if isnumeric(s.peek)
  fhan = s.get;
end
prob = s.get;
if ischar(s.peek) && strcmpi(s.peek, 'all')
  all = true;
end

if isfield(prob, 'signals')
  t        = coco_opts_tree();
  sigs     = prob.signals;
  signames = fieldnames(sigs);
  for i=1:numel(signames)
    signame = signames{i};
    path    = sigs.(signame).name;
    if all || has_slot(prob, signame)
      t = t.prop_set(path, 'owner', sigs.(signame).owner);
    end
    t = add_slots(prob, path, signame, t);
  end
  t.print_tree(fhan, 'signals', true);
end

end

function flag = has_slot(prob, signame)
flag = isfield(prob, 'slots') && isfield(prob.slots, signame);
end

function t = add_slots(prob, path, signame, t)
if has_slot(prob, signame)
  funcs = prob.slots.funcs;
  idx   = prob.slots.(signame);
  slots = {};
  for i=idx
    slotid = sprintf('%s@%s', funcs(i).identifyer, func2str(funcs(i).F));
    slots = [ slots slotid ]; %#ok<AGROW>
  end
  t = t.prop_set(path, 'slots', slots);
end
end
