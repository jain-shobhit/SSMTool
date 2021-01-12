function opts = state_BP_locate(opts)
%STATE_LOCATE_BOUNDARY_EVENTS  Locate events in current simplex.
%
%   OPTS = STATE_LOCATE_BOUNDARY_EVENTS(OPTS) assigns 'sol' to 'ptlist' and
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

BP_idx = opts.efunc.ev.BP_idx;
if isempty(BP_idx)
  opts.cont.state = 'MX_check';
else
  [opts opts.cseg] = opts.cseg.evlist_init(opts, BP_idx);
  if isempty(opts.cseg.evlist)
    opts.cont.state             = 'MX_check';
  else
    % opts.cseg.src_chart         = opts.cseg.ptlist{end};
    opts.cont.locate_next_state = 'MX_check';
    opts.cont.locate_warn_state = 'BP_warning';
    opts.cont.locate_add_state  = 'BP_add';
    opts.cont.state             = 'ev_locate';
  end
end
