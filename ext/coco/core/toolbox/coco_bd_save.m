function opts = coco_bd_save(opts)
  [opts fids data] = coco_emit(opts, 'save_bd');
  bd_data = [ fids data ]; %#ok<NASGU>
  version = 1; %#ok<NASGU>
	save(opts.run.bdfname, 'version', 'bd_data');
end
