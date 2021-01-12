function [prob, atlas, cseg, flush] = add_chart(atlas, prob, cseg)
if atlas.first
  chart = cseg.ptlist{1};
  cseg.curr_chart.ignore_evs = chart.ignore_evs;
  cseg.curr_chart.pt_type = chart.pt_type;
  cseg.curr_chart.ep_flag = chart.ep_flag;
  [prob, cseg]  = cseg.add_chart(prob, cseg.curr_chart);
  atlas.chart_list{1} = cseg.ptlist{2};
  flush = false;
  return
end
chart    = cseg.curr_chart;
chart.pt = chart.pt + 1;
if chart.pt >= atlas.PtMX(1)
  chart.pt_type = 'EP';
  chart.ep_flag = 1;
end
[prob, cseg]  = cseg.add_chart(prob, chart);
flush        = true;

cont = atlas.cont;

% check angle condition and refine step size if necessary
% This test is a modification of `arc_beta > prob.cont.arc_alpha' and
% allows for accepting charts with angles somewhat larger than arc_alpha.
% This modification significantly reduces the amount of unnecessarily
% rejected steps. The 'abs' is necessary to avoid complex numbers
% that occur sometimes due to roundoff errors.
chart1   = cseg.ptlist{1};
chart2   = cseg.ptlist{2};
t1       = chart1.t;
R1       = chart1.R;
arc_beta = abs(acos(t1' * chart2.t));
coco_log(prob, 2, cont.LogLevel, '%s: angle(t1,t0) = %.2e[DEG]\n', ...
  mfilename, 180 * arc_beta / pi);
if (chart2.pt>0) && (arc_beta > cont.h_fac_max * cont.arc_alpha)
  
  % reduce step size if possible and repeat continuation step
  if R1 > cont.h_min
    coco_warn(prob, 3, cont.LogLevel, ...
      'atlas: beta [%.4e] > h_fac_max * al_max, refining step size\n', ...
      180 * arc_beta / pi);
    
    chart1.R            = max(cont.h_min, cont.h_fac_min*R1);
    atlas.chart_list{1} = chart1;
    atlas.first         = true;
    flush               = false;
    return
    
  else % R1 <= cont.h_min
    coco_warn(prob, 3, cont.LogLevel, ...
      'atlas: minimum stepsize reached, but beta [%.4e] > h_fac_max * al_max\n', ...
      180 * arc_beta / pi);
  end
  
end

% copmpute new radius
if chart2.pt>0
  % This if-statement takes care of the case arc_beta==0, which produced an
  % error in some versions of Matlab.
  if cont.h_fac_max^2 * arc_beta < cont.arc_alpha
    h_fac = cont.h_fac_max;
  else
    h_fac = cont.arc_alpha / (sqrt(cont.h_fac_max) * arc_beta);
    h_fac = max(cont.h_fac_min, min(cont.h_fac_max, h_fac));
  end
  chart2.R = cont.ga * h_fac * R1;
  chart2.R = max( min(chart2.R, cont.h_max), cont.h_min);
end

% check residuum at new vertex
prob = coco_emit(prob, 'set_mode', 'check_res');
while chart2.R>cont.h_min
  x = chart2.x + chart2.TS*(chart2.R*chart2.s);
  [prob, chart2, f] = prob.efunc.F(prob, chart2, x);
  % coco_plot_F(prob, f); coco_plot_u(prob, x); coco_plot_u(prob, chart2.TS);
  if norm(f) > cont.MaxRes
    coco_warn(prob, 3, cont.LogLevel, ...
      'atlas: norm(f)=%.4e, refining step size\n', norm(f));
    h_fac   = max(cont.h_fac_min, cont.MaxRes/norm(f));
    chart2.R = max(cont.h_min, h_fac*chart2.R);
  else
    break;
  end
end

% copy chart back as end-point of ptlist
cseg.ptlist{2} = chart2;
end
