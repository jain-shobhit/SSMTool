function opts = state_ev_locate(opts)

events = opts.cont.events;
cseg   = opts.cseg;

%% check if we are in a reverse communication loop

if isfield(events, 'hanidx')
	
	events.call_count = events.call_count + 1;

	if events.call_count > events.callMX
		fprintf(2, '%s: warning: %s ''@%s'' %s [%d]\n%s: %s\n', ...
			mfilename, 'number of calls to event handler', ...
			func2str(opts.efunc.events(events.hanidx).han), ...
      'exceeded callMX', events.callMX, ...
			mfilename, 'locating events aborted');
		opts.cont.events            = rmfield(events, 'hanidx');
		opts.cont.locate_warn_state = events.warn_state;
		opts.cont.locate_add_state  = events.add_state;
		opts.cont.state             = 'ev_locate';
		return
	end
	
	switch events.state
	
		case 'init'
      evhan = opts.efunc.events(events.hanidx).han;
      data  = opts.efunc.events(events.hanidx).data;
      cseg  = opts.cseg;
			[data cseg events.msg] = evhan(opts, data, cseg, 'init', events.msg);
      opts.cseg = cseg;
      opts.efunc.events(events.hanidx).data = data;
			
			if events.call_count==1 && isfield(events.msg, 'callMX')
				events.callMX = events.msg.callMX;
			end
			
			switch events.msg.action

				case 'locate'
					events.evidx                = events.evgroup(events.msg.idx);
					events.point_type           = events.msg.point_type;
					events.state                = 'check';
					opts.cont.locate_warn_state = events.warn_state;
					opts.cont.locate_add_state  = 'ev_locate';
					
				case 'finish'
					opts.cont.events            = rmfield(events, 'hanidx');
					opts.cont.locate_warn_state = events.warn_state;
					opts.cont.locate_add_state  = events.add_state;
					opts.cont.state             = 'ev_locate';
					return
          
				case 'warn'
					opts.cont.events            = events;
					opts.cont.locate_warn_state = events.warn_state;
					opts.cont.locate_add_state  = events.add_state;
					opts.cont.state             = events.warn_state;
					return
          
        otherwise
          error('%s: illegal action requested by event handler ''%s''', ...
            mfilename, func2str(evhan));
			end
		
		case 'check'
      evhan = opts.efunc.events(events.hanidx).han;
      data  = opts.efunc.events(events.hanidx).data;
      cseg  = opts.cseg;
			[data cseg events.msg] = evhan(opts, data, cseg, 'check', events.msg);
      opts.cseg = cseg;
      opts.efunc.events(events.hanidx).data = data;

			switch events.msg.action
				
				case 'add'
					events.state                = 'init';
          events.point_type           = events.msg.point_type;
					opts.cont.events            = events;
					opts.cont.locate_warn_state = events.warn_state;
					opts.cont.locate_add_state  = events.add_state;
					opts.cont.state             = events.add_state;
					return
					
				case 'reject'
					events.state                = 'init';
					opts.cont.events            = events;
					opts.cont.locate_warn_state = events.warn_state;
					opts.cont.locate_add_state  = events.add_state;
					opts.cont.state             = 'ev_locate';
					return
          
        otherwise
          error('%s: illegal action requested by event handler ''%s''', ...
            mfilename, func2str(evhan));
			end
	end
	
%% initialise new event group

else

%% leave event handling if all events located

	if isempty(cseg.evlist)
		opts.cont.state = opts.cont.locate_next_state;
		return
	end

%% compute indices of hit event group

	evlist  = cseg.evlist;
	evidx   = evlist(1);
	evgroup = opts.efunc.ev.evgroup{evidx};
  hanidx  = opts.efunc.ev.idx(evidx);

	if ~isempty(evgroup)
		evidx = evgroup;
	end

%% remove current event group from event list

	for i=1:numel(evidx)
		evlist(evlist==evidx(i)) = [];
	end

	events.evidx = evidx;
	cseg.evlist  = evlist;
	
	if ~isempty(opts.efunc.events(hanidx).han)
		msg.u0         = cseg.u0;
		msg.u1         = cseg.u1;
		msg.e0         = cseg.e0(evidx);
		msg.e1         = cseg.e1(evidx);
		msg.pidx       = opts.efunc.ev.pidx(evidx);
		msg.midx       = opts.efunc.ev.midx(evidx);
		msg.p0         = cseg.p0(msg.midx);
		msg.p1         = cseg.p1(msg.midx);
		msg.pars       = coco_idx2par(opts, msg.pidx);
		msg.evidx      = evidx;
    msg.action     = '';
    msg.idx        = [];
    msg.point_type = '??';
		
    events.hanidx     = hanidx;
		events.msg        = msg;
		events.evgroup    = evidx;
		events.warn_state = opts.cont.locate_warn_state;
		events.add_state  = opts.cont.locate_add_state;
		events.call_count = 0;
		events.callMX     = 100;

		events.state     = 'init';
		opts.cont.state  = 'ev_locate';
		opts.cont.events = events;
    opts.cseg        = cseg;
		return
	end

%% end of check if we are in a reverse communication loop

end

%% initialise cache used by subdivision methods

cseg.ev0  = cseg.e0(events.evidx);
cseg.v0   = cseg.u0;
cseg.h0   = 0;
cseg.ev1  = cseg.e1(events.evidx);
cseg.v1   = cseg.u1;
cseg.h1   = 1;
opts.cseg = cseg;

%% dispatch computation of event to appropriate algorithm

if numel(events.evidx)>1
	opts.cont.state = 'ev_locate_multi';
else
	switch opts.efunc.ev.par_type{events.evidx}

		case 'continuation'
			opts.cont.state = 'ev_locate_cont';

		case 'regular'

			% initialise subdivision method, placing this initialisation here
			% implies that ev_locate_reg is always entered after correction
			e0 = cseg.ev0;
			e1 = cseg.ev1;

			v0 = cseg.v0;
			h0 = cseg.h0;
			v1 = cseg.v1;
			h1 = cseg.h1;

			la1 = cseg.la1;
			la2 = cseg.la2;

      if e0*e1==0
        % go straight to add_<PointType>
        if e0==0
          cseg.curr_chart = cseg.ptlist{1};
        else
          cseg.curr_chart = cseg.ptlist{end};
        end
        opts.cont.state = opts.cont.locate_add_state;
      else
        % compute first subdivision point
        if abs(e0)<=abs(e1)
          v = la2*v0+la1*v1;
          h = la2*h0+la1*h1;
        else
          v = la1*v0+la2*v1;
          h = la1*h0+la2*h1;
        end
        cseg.h = h;
        
        % interpolate
        [opts cseg chart hh] = cseg.chart_at(opts, h, events.evidx);
        
        % update projection condition
        opts = coco_emit(opts, 'update_h', hh);
        
        % initialise nwtn
        [opts chart accept x] = opts.corr.init(opts, chart, v);
        chart.x               = x;
        cseg.curr_chart       = chart;
        
        opts.cseg = cseg;
        
        % go to next state
        if accept
          % go to locate_reg
          opts.cont.state      = 'ev_locate_reg';
        else
          % next state is correct, then go to locate_reg
          opts.cont.state      = 'co_correct';
          opts.cont.next_state = 'ev_locate_reg';
          opts.cont.err_state  = opts.cont.locate_warn_state;
        end
      end
      
    case 'singular'
			opts.cont.state = 'ev_locate_sing';
	end
end

%% update opts.events
opts.cont.events = events;
