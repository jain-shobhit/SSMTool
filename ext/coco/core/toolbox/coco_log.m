function fhan = coco_log(opts, prio, LogLevel, varargin)

if nargin==3
  if isfield(opts, 'run') && prio<=max(opts.run.logPrioMN, LogLevel(1))
    fhan = opts.run.loghan;
  else
    fhan = [];
  end
else
  if isfield(opts, 'run')
    run      = opts.run;
    LogLevel = max(run.logPrioMN, LogLevel(1));
    if prio<=LogLevel
      fprintf(run.loghan, varargin{:});
    end
  end
end

end
