function [opts, all_funcs] = coco_save_funcs(opts)
all_funcs.efunc   = opts.efunc;
if isfield(opts, 'adjoint')
  all_funcs.adjoint = opts.adjoint;
end
all_funcs.slots   = opts.slots;
if isfield(opts.efunc, 'events')
  all_funcs.events = opts.efunc.events;
  all_funcs.ev     = opts.efunc.ev;
end
all_funcs.ptrs = coco_func_data.pointers('copy');

% call toolbox data copy functions

nefuncs = numel(opts.efunc.funcs);
nsfuncs = numel(opts.slots.funcs);
if isfield(opts.efunc, 'events')
  nhfuncs = numel(opts.efunc.events);
else
  nhfuncs = 0;
end

for i=1:nefuncs
  data = opts.efunc.funcs(i).data;
  fid  = opts.efunc.funcs(i).identifyer; %#ok<NASGU> % for debugging only
  copy = opts.efunc.funcs(i).copy;
  if ~isempty(copy)
    opts.efunc.funcs(i).data = copy(opts, data);
  end
end

for i=1:nsfuncs
  data = opts.slots.funcs(i).data;
  fid  = opts.slots.funcs(i).identifyer; %#ok<NASGU> % for debugging only
  copy = opts.slots.funcs(i).copy;
  if ~isempty(copy)
    opts.slots.funcs(i).data = copy(opts, data);
  end
end

for i=1:nhfuncs
  data = opts.efunc.events(i).data;
  fid  = opts.efunc.events(i).name{1}; %#ok<NASGU> % for debugging only
  copy = opts.efunc.events(i).copy;
  if ~isempty(copy)
    opts.efunc.events(i).evdata{i} = copy(opts, data);
  end
end

end
