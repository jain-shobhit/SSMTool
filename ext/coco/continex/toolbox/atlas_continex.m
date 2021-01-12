classdef atlas_continex < AtlasBase
  % atlas_continex  Basic atlas code with non-linear arc-length condition.
  
  properties (Access=private)
    base_chart = struct();
    cont       = struct();
    pts        = [];
    gapnum     = 0;
  end
  
  methods (Access=private)
    function atlas = atlas_continex(prob, dim)
      assert(dim==1, '%s: wrong manifold dimension', mfilename);
      atlas      = atlas@AtlasBase(prob);
    end
  end
  
  methods (Static)
    function [prob cont atlas] = create(prob, cont, dim)
      prob       = coco_add_signal(prob, 'continex_force_FDM', mfilename);
      atlas      = atlas_continex(prob, dim);
      atlas.cont = atlas.get_settings(cont);
      atlas.cont = cseg_continex.get_settings(atlas.cont);
      cont       = atlas.cont;
      cont.corr  = struct(); % we take care of post init settings
      prob       = cseg_continex.add_prcond(prob, cont, dim);
      
      func = @(x) ( isfield(cont, x) && ~isempty(cont.(x)) );
      if all(cellfun(func, {'phan' 'xidx' 'yidx'}))
        prob = coco_add_slot(prob, 'atlas.plot_init', @plot_init, [], 'corr_begin');
        prob = coco_add_slot(prob, 'atlas.plot_cseg', @plot_cseg, [], 'corr_step');
        prob = coco_add_slot(prob, 'atlas.plot_smpl', @plot_smpl, [], 'corr_sample');
      end
    end
  end
  
  methods (Static, Access=private)
    function cont = get_settings(cont)
      defaults.PtMX  = 50   ; % number of continuation steps
      defaults.h_min = 0.05 ; % minimum step size
      defaults.h0    = 0.05 ; % initial step size
      defaults.h_max = 0.1  ; % maximum step size
      defaults.NPtAv = 1    ; % number of points for fitting TS
      defaults.GAPMX = inf  ; % terminate after GAPMX consecutive GAP points
      defaults.TrMX  = 1    ; % number of trials
      
      defaults.phan = []; % debug plot axis handle
      defaults.xidx = []; % index of x-axis
      defaults.yidx = []; % index of y-axis
      
      defaults.corrector = 'continex';
      
      cont = coco_merge(defaults, cont);
      
      % bug: remove assertion after fixing run_cont etc.
      assert(isscalar(cont.TrMX) && cont.TrMX>=1);
      
      if isfield(cont, 'ItMX')
        cont.PtMX = cont.ItMX;
      end
    end
  end
  
  methods (Access=public)
    
    function [prob atlas cseg correct] = init_prcond(atlas, prob, chart)
      atlas.cont.corr = init_trial_opts(prob, atlas);
      chart.R       = 0;
      chart.pt      = -1;
      chart.pt_type = 'IP';
      chart.ep_flag = 1;
      [prob cseg]   = cseg_continex.create_initial(prob, chart);
      correct       = cseg.correct;
    end
    
    function [prob atlas cseg flush] = init_atlas(atlas, prob, cseg)
      chart           = cseg.curr_chart;
      chart.R         = atlas.cont.h0;
      chart.s         = sign(atlas.cont.PtMX);
      atlas.cont.PtMX = abs(atlas.cont.PtMX);
      chart.pt        = 0;
      chart.pt_type   = 'EP';
      chart.ep_flag   = 1;
      atlas.pts       = chart.x;
      [prob cseg]     = cseg.add_chart(prob, chart, atlas.pts);
      flush           = true;
    end
    
    function [prob atlas cseg] = flush(atlas, prob, cseg)
      [prob atlas cseg] = atlas.flush@AtlasBase(prob, cseg);
      if cseg.Status == cseg.CurveSegmentOK
        atlas.base_chart = cseg.ptlist{end};
        if atlas.base_chart.pt>=atlas.cont.PtMX
          cseg.Status = cseg.BoundaryPoint;
        end
      end
    end
    
    function [prob atlas cseg correct] = predict(atlas, prob, cseg)
      correct = true;
      if isempty(cseg)
        chart     = atlas.base_chart;
        chart.gap = false;
        if chart.pt>0
          prcond = struct('x', chart.x, 'TS', chart.TS, ...
            's', chart.s, 'h', atlas.cont.h_max);
          corr        = atlas.cont.corr{1};
          prob        = coco_emit(prob, 'shape_set_scale', corr.scale);
          [prob cseg] = cseg_continex.create(prob, chart, prcond);
          cseg.corr   = corr;
        else
          prcond = struct('x', chart.x, 'TS', chart.TS, ...
            's', chart.s, 'h', atlas.cont.h0);
          [prob cseg] = cseg_continex.create(prob, chart, prcond);
        end
      else
        prob = coco_emit(prob, 'shape_set_scale', cseg.corr.scale);
        [prob cseg] = cseg.reset(prob);
      end
      if cseg.base_chart.pt>0
        cseg.corr.MaxAbsDist = cseg.MaxAbsDist;
        cseg.corr.MaxAbsStep = cseg.MaxAbsDist/cseg.corr.MADSteps;
        prob.corr = prob.corr.set_opts(prob.corr, cseg.corr);
        cseg.curr_chart.gap = cseg.corr.GAP;
        if cseg.corr.Jac
          prob = coco_emit(prob, 'continex_force_FDM');
        end
      end
    end
    
    function [prob atlas cseg flush] = add_chart(atlas, prob, cseg)
      chart    = cseg.curr_chart;
      chart.pt = chart.pt + 1;
      chart.R  = norm(chart.x-cseg.base_chart.x);
      if chart.gap
        chart.pt_type = 'GAP';
        atlas.gapnum  = atlas.gapnum+1;
        if chart.pt >= atlas.cont.PtMX
          [prob cseg] = cseg.add_chart(prob, chart, atlas.pts);
        end
      else
        atlas.gapnum = 0;
        atlas.pts(:,end+1) = chart.x;
        % bug: reorder points wrt. average direction
        if size(atlas.pts,2)>atlas.cont.NPtAv
          atlas.pts(:,1) = [];
        end
      end
      if (atlas.cont.GAPMX>0) && (atlas.gapnum>=atlas.cont.GAPMX)
        chart.pt_type = 'GMX';
        chart.ep_flag = 2;
        cseg.Status   = cseg.CurveSegmentCorrupted;
      elseif chart.pt >= atlas.cont.PtMX
        chart.pt_type = 'EP';
        chart.ep_flag = 1;
      end
      [prob cseg] = cseg.add_chart(prob, chart, atlas.pts);
      flush       = true;
    end
    
    function [prob atlas cseg predict] = refine_step(atlas, prob, cseg)
      predict = false;
      if cseg.base_chart.pt>0
        cseg.trial = cseg.trial+1;
        if cseg.trial<=atlas.cont.TrMX
          cseg.corr = atlas.cont.corr{cseg.trial};
          predict   = true;
          coco_warn(prob, 3, atlas.cont.LogLevel, ...
            'no convergence, trying again, trial %d of %d\n', ...
            cseg.trial, atlas.cont.TrMX);
        end
      end
    end
    
  end
  
end

function corr = init_trial_opts(prob, atlas)
cont     = atlas.cont;
scales   = linspace(1, cont.h_min/cont.h_max, cont.TrMX);
defaults = struct('Jac', false, 'GAP', false, 'scale', 1);
defaults = coco_merge(defaults, prob.corr.get_opts(prob.corr));
corr     = cell(1, cont.TrMX);
corr{1}  = coco_merge(defaults, cont.corr);
if isfield(cont, 'corr1')
  corr{1} = coco_merge(corr{1}, cont.('corr1'));
  if isfield(cont.('corr1'), 'scale')
    scales(1:cont.TrMX) = cont.('corr1').scale;
  end
end
for i=2:cont.TrMX
  corr{i} = corr{i-1};
  corr{i}.scale = scales(i);
  corr_i = sprintf('corr%d', i);
  if isfield(cont, corr_i)
    corr{i} = coco_merge(corr{i}, cont.(corr_i));
    if isfield(cont.(corr_i), 'scale')
      scales(i:cont.TrMX) = cont.(corr_i).scale;
    end
  end
end
end

function data = plot_init(prob, data, corrector, chart, x0, varargin) %#ok<INUSL>
atlas  = prob.atlas;
cont   = atlas.cont;
phan   = cont.phan;
xidx   = cont.xidx;
yidx   = cont.yidx;
cseg   = prob.cseg;
prcond = cseg.prcond;

cla(prob.atlas.cont.phan, 'reset');
hold(phan, 'on');

u = atlas.pts;
if ~isempty(u)
  x = u(xidx,:);
  y = u(yidx,:);
  plot(phan, x, y, 'g*');
  
  uu = repmat(u(:,end), 1, size(u,2));
  t  = prcond.TS*prcond.s;
  v  = uu + kron(t, t'*(u-uu));
  x  = v(xidx,:);
  y  = v(yidx,:);
  plot(phan, x, y, 'ko');
  v  = [v(:,1)-prcond.h*t v v(:,end)+prcond.h*t];
  x  = v(xidx,:);
  y  = v(yidx,:);
  plot(phan, x, y, 'k-');
end

if chart.pt>0
  cseg.plot_prcond(prob, phan, xidx, yidx);
  mrk = prob.corr.print.mrk;
  sample = strcmp(mrk, regexp(mrk, '\d+', 'match', 'once'));
  if sample
    plot(phan, x0(xidx), x0(yidx), 'kp', x0(xidx), x0(yidx), 'ro');
  else
    plot(phan, x0(xidx), x0(yidx), 'ro');
  end
  prob.cseg.plot_point(prob, phan, x0, xidx, yidx, sample);
end

grid(phan, 'on');
hold(phan, 'off');
axis(phan, 'equal');

end

function [data stop msg] = plot_cseg(prob, data, corrector, chart, x, varargin) %#ok<INUSL>
[stop msg] = deal(false, '');

atlas = prob.atlas;
cont  = atlas.cont;
phan  = cont.phan;
xidx  = cont.xidx;
yidx  = cont.yidx;

hold(phan, 'on');
if chart.pt>0
  mrk = prob.corr.print.mrk;
  sample = strcmp(mrk, regexp(mrk, '\d+', 'match', 'once'));
  if sample
    plot(phan, x(xidx), x(yidx), 'kp', x(xidx), x(yidx), 'ro');
  else
    plot(phan, x(xidx), x(yidx), 'ro');
  end
  prob.cseg.plot_point(prob, phan, x, xidx, yidx, sample);
end
hold(phan, 'off');

end

function data = plot_smpl(prob, data, x, varargin) %#ok<INUSD>

atlas = prob.atlas;
cont  = atlas.cont;
phan  = cont.phan;
xidx  = cont.xidx;
yidx  = cont.yidx;
cseg  = prob.cseg;

hold(phan, 'on');
if cseg.curr_chart.pt>0
  cseg.plot_prcond(prob, phan, xidx, yidx);
end
hold(phan, 'off');

end
