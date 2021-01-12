classdef EventLocator
  
  properties % public properties
    
    % properties for event handling, these cache special properties of the
    % end points of a curve segment for faster access
    u0, p0, e0, ign0, u1, p1, e1, ign1;

    % properties for event handling, these cache special properties
    % computed during event handling
    ev0, ev1, evlist;
    v0, v1, h0, h1, h;
    
  end
  
  properties (Constant = true) % predefined constants
    gm  = (1+sqrt(5))/2;
    la1 = 1-2/(1+sqrt(5));
    la2 = 2/(1+sqrt(5));
  end
  
  methods (Abstract=true)
    [opts cseg chart h] = chart_at     (cseg, opts, la, chart, evidx)
    [opts cseg idx]     = insert_chart (cseg, opts, chart)
    [opts cseg]         = eval_p       (cseg, opts, evidx)
  end
  
  methods % interface methods to event states of FSM
    
    function [opts cseg] = ev_init(cseg, opts)
      
      % initialise event fields
      for i=1:numel(cseg.ptlist)
        chart          = cseg.ptlist{i};
        [opts chart.e] = opts.efunc.events_F(opts, chart.p);
        cseg.ptlist{i} = chart;
      end
      
      % initialise cseg for event handling
      chart0  = cseg.ptlist{1};
      cseg.u0 = chart0.x;
      cseg.p0 = chart0.p;
      cseg.e0 = chart0.e;
      if isfield(chart0, 'ignore_evs')
        cseg.ign0 = chart0.ignore_evs;
      else
        cseg.ign0 = [];
        cseg.ptlist{1}.ignore_evs = [];
      end
      if ~isfield(chart0, 'ignore_at')
        cseg.ptlist{1}.ignore_at  = false;
      end
      
      chart1  = cseg.ptlist{end};
      cseg.u1 = chart1.x;
      cseg.p1 = chart1.p;
      cseg.e1 = chart1.e;
      if isfield(chart1, 'ignore_evs')
        cseg.ign1 = chart1.ignore_evs;
      else
        cseg.ign1 = [];
        cseg.ptlist{end}.ignore_evs = [];
      end
      if ~isfield(chart1, 'ignore_at')
        cseg.ptlist{end}.ignore_at  = false;
      end
      
    end
    
    function [opts cseg] = evlist_init(cseg, opts, idx)
      
      % initialise event list
      cseg.ev0 = cseg.e0;
      cseg.ev1 = cseg.e1;
      
      evmask         = false(numel(cseg.ev1),1);
      evmask(idx)    = true;
      
      ignore         = intersect(idx, cseg.ign0);
      evmask(ignore) = false;
      if cseg.ptlist{  1}.ignore_at
        ignore         = idx(cseg.ev0(idx) == 0);
        evmask(ignore) = false;
      end
      
      ignore           = intersect(idx, cseg.ign1);
      evmask(ignore)   = false;
      if cseg.ptlist{end}.ignore_at
        ignore         = idx(events.ev1(idx) == 0);
        evmask(ignore) = false;
      end
      
      evcrossed   = cseg.ev0.*cseg.ev1<=0;
      % bug: include log message about constant but crossed events
      evconst     = abs(cseg.ev0-cseg.ev1)<=10*opts.cont.TOL;
      cseg.evlist = find(evmask & evcrossed & ~evconst);
      cseg.evlist = cseg.evlist(:)';
      
    end
    
    function [opts cseg] = add_BP(cseg, opts, chart)
      
      [opts cseg idx] = cseg.add_SP(opts, chart);
      
      cseg.ptlist = cseg.ptlist(1:idx);
      
      chart1  = cseg.ptlist{end};
      cseg.u1 = chart1.x;
      cseg.p1 = chart1.p;
      cseg.e1 = chart1.e;
      if isfield(chart1, 'ignore_evs')
        cseg.ign1 = chart1.ignore_evs;
      else
        cseg.ign1 = [];
        cseg.ptlist{end}.ignore_evs = [];
      end
      if ~isfield(chart1, 'ignore_at')
        cseg.ptlist{end}.ignore_at  = false;
      end
      
    end
    
    function [opts cseg idx] = add_SP(cseg, opts, chart)
      [opts cseg idx]  = cseg.insert_chart(opts, chart);
      chart            = cseg.ptlist{idx};
      %[opts chart]     = opts.cseg.update_det(opts, chart); % Harry added: Compute determinant of restricted problem
      [opts chart.e]   = opts.efunc.events_F(opts, chart.p);
      cseg.ptlist{idx} = chart;
    end
    
  end
  
end
