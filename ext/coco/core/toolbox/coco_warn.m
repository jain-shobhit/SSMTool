function coco_warn(opts, prio, LogLevel, varargin)

if isstruct(opts) && isfield(opts, 'run')
  run      = opts.run;
  LogLevel = max(run.logPrioMN, LogLevel(1));
  fprintf(run.loghan, '\nwarning: ');
  fprintf(run.loghan, varargin{:});
end

if prio<=LogLevel(1)
  fprintf(2, '\nwarning: ');
  fprintf(2, varargin{:});
end

end
