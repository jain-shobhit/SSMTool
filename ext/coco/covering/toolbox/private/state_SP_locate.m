function opts = state_SP_locate(opts)
%STATE_COMPUTE_EVENTS  Locate events in current simplex.
%
%   OPTS = STATE_COMPUTE_EVENTS(OPTS) assigns 'sol' to 'ptlist' and
%   instantiates an empty 'spt' property if necessary. If the maximum
%   number of iterations along the solution branch in the current direction
%   has been reached, then this property is set to 'EP'. (Call to ptlist(1)
%   seems unnecessary). The state is subsequently set to HANDLE_EVENTS.
%
%   To do: This is confusing. Why is ptlist erased here. It was just given
%   a value in check_solution in case the boundary was reached. On the
%   other hand, if the boundary was not reached, ptlist doesn't even exist.
%
%   See also:
%

% compute list of events

SP_idx = opts.efunc.ev.SP_idx;
if isempty(SP_idx)
  opts.cont.state = 'co_flush';
else
  [opts opts.cseg] = opts.cseg.evlist_init(opts, SP_idx);
  if isempty(opts.cseg.evlist)
    opts.cont.state             = 'co_flush';
  else
    % opts.cseg.src_chart         = opts.cseg.ptlist{end};
    opts.cont.locate_next_state = 'co_flush';
    opts.cont.locate_warn_state = 'SP_warning';
    opts.cont.locate_add_state  = 'SP_add';
    opts.cont.state             = 'ev_locate';
  end
end
