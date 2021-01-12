function opts = coco_add_pending(opts, satisfied)

if nargin<2
  satisfied = {};
end

if ischar(satisfied)
  satisfied = { satisfied };
end

if ~isfield(opts, 'efunc')
  opts.efunc = efunc_new([]);
end
if ~isfield(opts.efunc, 'identifyers')
  opts.efunc = efunc_new(opts.efunc);
end

% disable calling add_pending while resolving dependencies
opts.efunc.add_pending = false;

% add pending functions
func_added = true;
while(func_added)
  func_added = false;
  for func = 1:size(opts.efunc.pending,1)
    deplist        = opts.efunc.pending{func,1};
    fids           = [opts.efunc.identifyers satisfied];
    deps_satisfied = true;
    for i = 1:numel(deplist)
      if ~any(strcmp(deplist{i}, fids))
        deps_satisfied = false;
        break
      end
    end
    if deps_satisfied
      ctor = opts.efunc.pending{func,2};
      args = opts.efunc.pending{func,3};
      opts = ctor(opts, args{:});
      opts.efunc.pending(func,:) = [];
      func_added = true;
      break;
    end
  end
end

% enable calling add_pending
opts.efunc.add_pending = true;

end
