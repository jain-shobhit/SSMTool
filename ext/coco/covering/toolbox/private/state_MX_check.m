function opts = state_MX_check(opts)
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

% Check for terminal events

MX_idx = opts.efunc.ev.MX_idx;
if isempty(MX_idx)
  opts.cont.state = 'SP_locate';
else
  [opts opts.cseg] = opts.cseg.evlist_init(opts, MX_idx);
  if isempty(opts.cseg.evlist)
    opts.cont.state = 'SP_locate';
  else
    evlist = opts.cseg.evlist;
    opts.cseg.ptlist{end}.pt_type = opts.efunc.ev.point_type{evlist(1)};

    % mark curve segment as invalid
    opts.cseg.Status = opts.cseg.TerminalEventDetected;
    
    % next state is flush
    opts.cont.state = 'co_flush';
  end
end
