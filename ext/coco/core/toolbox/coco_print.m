function coco_print(opts, prio, varargin)

if isfield(opts, 'run')
  run = opts.run;
  if prio<=run.scrPrioMN
    fprintf(run.scrhan, varargin{:});
    fprintf(1, varargin{:});
  end
  fprintf(run.loghan, varargin{:});
else
  fprintf(1, varargin{:});
end

end
