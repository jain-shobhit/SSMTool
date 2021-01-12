function [prob, cont, corr] = corr_fsolve(prob, cont, func)

% set up Matlab's fsolve
[prob, corr] = get_settings(prob);

corr.FDF      = prob.(func).FDF;
corr.solve    = @solve;
corr.init     = @init;
corr.step     = @step;
corr.set_opts = @set_opts;
corr.get_opts = @get_opts;

prob = coco_add_signal(prob, 'corr_begin', mfilename);
prob = coco_add_signal(prob, 'corr_step', mfilename);
prob = coco_add_signal(prob, 'corr_end', mfilename);

prob = coco_add_signal(prob, 'corr_print', mfilename);

end

function [prob, corr] = get_settings(prob)

defaults.ItMX = 10      ; % max. number of iterations
defaults.opts = optimset; % options passed to fsolve

corr = coco_get(prob, 'corr');
corr = coco_merge(defaults, corr);

defaults.ResTOL = corr.TOL; % condition: ||f(x)||<=ResTOL

corr = coco_merge(defaults, corr);

% copy options for fsolve
if isempty(optimget(corr.opts, 'Jacobian'))
  corr.opts = optimset(corr.opts, 'Jacobian', 'on');
end
if isempty(optimget(corr.opts, 'MaxIter'))
  corr.opts = optimset(corr.opts, 'MaxIter', corr.ItMX);
end
if isempty(optimget(corr.opts, 'TolX'))
  corr.opts = optimset(corr.opts, 'TolX', corr.TOL);
end
if isempty(optimget(corr.opts, 'TolFun'))
  corr.opts = optimset(corr.opts, 'TolFun', corr.ResTOL);
end

corr.filter = {'LogLevel' 'TOL' 'ResTOL' 'ItMX' 'lsol' 'opts'};

end

function corr = set_opts(corr, settings)
corr = coco_merge(corr, settings, corr.filter);
end

function settings = get_opts(corr)
settings = coco_merge(struct(), corr, corr.filter);
end

function [prob, chart, x] = solve(prob, chart, x0)

[prob, chart, accept, x] = init(prob, chart, x0);
while ~accept
	[prob, chart, accept, x] = step(prob, chart);
end

end

function [prob, chart, accept, x0] = init(prob, chart, x0)

corr = prob.corr;

if 3 <= prob.run.scrPrioMN
  corr.opts = optimset(corr.opts, 'Diagnostics', 'on'   );
  corr.opts = optimset(corr.opts, 'Display'    , 'iter' );
elseif 2 <= prob.run.scrPrioMN
  corr.opts = optimset(corr.opts, 'Diagnostics', 'off'  );
  corr.opts = optimset(corr.opts, 'Display'    , 'iter' );
else
  corr.opts = optimset(corr.opts, 'Diagnostics', 'off'  );
  corr.opts = optimset(corr.opts, 'Display'    , 'off'  );
end

corr.x = x0;
[prob, chart, f] = corr.FDF(prob, chart, x0);
accept = (norm(f)<corr.ResTOL);

prob.corr = corr;

prob = coco_emit(prob, 'corr_begin', 'fsolve', chart, x0);
if accept
  prob = coco_emit(prob, 'corr_end', 'fsolve', 'accept', chart, x0);
end

end

function [prob, chart, accept, x] = step(prob, chart)

corr = prob.corr;

corr.opts   = optimset(corr.opts, 'OutputFcn', @OFunc);
emsg        = ''; % set in call output function below
[x, fval, exitflag, output] = fsolve(@func, corr.x, corr.opts); %#ok<ASGLU>
corr.x      = x;
corr.accept = ( exitflag ==  1 );

prob.corr = corr;
accept    = corr.accept;

if accept
  prob = coco_emit(prob, 'corr_end', 'fsolve', 'accept', chart, x);
elseif exitflag == -1
  prob = coco_emit(prob, 'corr_end', 'fsolve', 'stop');
	errmsg.identifier = 'CORR:Stop';
	errmsg.message = sprintf(...
    '%s: stop requested by slot function(s)\n%s', ...
    mfilename, emsg);
	errmsg.FID = 'FSOLVE';
	errmsg.ID  = 'MX';
	error(errmsg);
else
  prob = coco_emit(prob, 'corr_end', 'fsolve', 'fail');
	errmsg.identifier = 'CORR:NoConvergence';
	errmsg.message = output.message;
	errmsg.FID = 'FSOLVE';
	errmsg.ID  = 'MX';
	error(errmsg);
end

  function varargout = func(x)
    [prob, chart, varargout{1:nargout}] = corr.FDF(prob, chart, x);
  end

  function stop = OFunc(x, optimValues, state)
    corr = prob.corr;
    
    stop   = false;
    corr.x = x;
    
    switch state
      case 'init'
        prob = coco_emit(prob, 'corr_begin', 'fsolve', ...
          x, optimValues, state);
      case 'interrupt'
        [prob, fids, stop, msg] = coco_emit(prob, 'corr_step', 'fsolve', ...
          x, state, optimValues, state);
        stop = cell2mat(stop);
        if ~isempty(stop) && any(stop)
          emsg = '';
          for idx=find(stop)
            emsg = sprintf('%s%s: %s\n', emsg, fids{idx}, msg{idx});
          end
          stop=true;
        else
          stop=false;
        end
    end
    
    prob.corr = corr;
  end

end
