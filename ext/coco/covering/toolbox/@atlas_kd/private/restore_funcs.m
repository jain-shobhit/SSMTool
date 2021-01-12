function opts = restore_funcs(opts, all_funcs)

nefuncs = numel(opts.efunc.funcs);
nsfuncs = numel(opts.slots.funcs);
if isfield(opts.efunc, 'events')
  nhfuncs = numel(opts.efunc.events);
else
  nhfuncs = 0;
end

for i=1:nefuncs
  if isa(opts.efunc.funcs(i).data, 'coco_func_data')
    opts.efunc.funcs(i).data.pr = all_funcs.efunc.funcs(i).data.pr;
    opts.efunc.funcs(i).data.sh = all_funcs.efunc.funcs(i).data.sh;
    all_funcs.efunc.funcs(i).data = opts.efunc.funcs(i).data;
  end
end

for i=1:nsfuncs
  if isa(opts.slots.funcs(i).data, 'coco_func_data')
    opts.slots.funcs(i).data.pr = all_funcs.slots.funcs(i).data.pr;
    opts.slots.funcs(i).data.sh = all_funcs.slots.funcs(i).data.sh;
    all_funcs.slots.funcs(i).data = opts.slots.funcs(i).data;
  end
end

for i=1:nhfuncs
  if isa(opts.efunc.events(i).data, 'coco_func_data')
    opts.efunc.events(i).data.pr = all_funcs.efunc.events(i).data.pr;
    opts.efunc.events(i).data.sh = all_funcs.efunc.events(i).data.sh;
    all_funcs.efunc.events(i).data = opts.efunc.events(i).data;
  end
end

opts.efunc   = all_funcs.efunc;
if isfield(all_funcs, 'adjoint')
  opts.adjoint = all_funcs.adjoint;
end
opts.slots   = all_funcs.slots;
if isfield(all_funcs, 'events')
  opts.efunc.events = all_funcs.events;
  opts.efunc.ev     = all_funcs.ev;
else
  if isfield(opts.efunc, 'events')
    opts.efunc = rmfield(opts.efunc, {'ev' 'events'});
  end
end

coco_func_data.pointers('restore', all_funcs.ptrs);

end
