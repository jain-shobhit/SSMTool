function prob = coco_add_event(varargin)
%COCO_ADD_EVENT   Add event to parameter.
%
% PROB = COCO_ADD_EVENT(PROB, EV_HAN, EV_SIG)
% Add events to a parameter. The arguments EV_HAN specify how to handle the
% events and the event signature EV_SIG defines the type of the event and
% associates event values to a monitor function.
%
% EV_HAN = { POINT_TYPE }
% Mark the event as type POINT_TYPE in the bifurcation diagram (column
% TYPE). The type identifyer POINT_TYPE is a string and should be of length
% 2-4. The same type identifyer can be used for several events.
%
% EV_HAN = { @FUNC DATA [@COPY] }
% Associate an event handler instance with the event. FUNC is a function
% handle to the handler, DATA is function data of the instance and COPY
% defines a function that creates a deep copy of DATA. COPY is usually not
% necessary. It should be provided if DATA is an instance of a handle class
% different from coco_func_data. The calling syntax of FUNC and COPY is
% explained below.
%
% EV_SIG = { [('SP'|'special point')] SP_SIG }
% SP_SIG = { PAR ['='] VALS } | GR_SIG
% GR_SIG = { ( PAR ['='] VAL )... }
% Event is a special point. Using GR_SIG defines a group event. For group
% events it is possible to associate parameters that do not trigger the
% event by assigning the value NaN to the parameter. This can be useful to
% make values of additional monitor functions available to an event
% handler.
%
% EV_SIG = { ('BP'|'boundary') BP_SIG }
% BP_SIG = { PAR ('<'|'>') VAL } | SP_SIG
% Event defines a computational boundary. The continuation will stop at the
% computational boundary. For '<' and '>' signatures, the event is
% triggered if PAR becomes less than or greater than VAL, respectively.
%
% EV_SIG = { ('MX'|'terminate') MX_SIG }
% MX_SIG = { PAR ['<'|'>'|'='] VAL }
% Event is a terminal event. The continuation will discard a curve segment
% containing a terminal event. For '<' and '>' signatures, the event is
% triggered if PAR becomes less than or greater than VAL, respectively.
%
% --------------------------------------------------
% [DATA CSEG MSG] = FUNC(PROB, DATA, CSEG, CMD, MSG)
% Call event handler FUNC following a reverse communication protocoll. CMD
% is either 'init' or 'check' and MSG is a structure that can be used
% freely by the event handler to record the history of calls for a
% particular special point to allow history-dependent execution. The
% following fields of the structure MSG are used by the default finite
% state mechine:
%
%   MSG.callMX     - Set maximum number of repeated calls to FUNC, this
%                    field is only used after the first call to FUNC for
%                    each point.
%   MSG.action     - Define action of event location after call to FUNC.
%                    Possible actions after CMD = 'init': 'locate',
%                    'finish', 'warn'. Possible actions after CMD =
%                    'check': 'add', 'reject'.
%   MSG.point_type - Type of special point. Can be modified by the event
%                    handler. MSG.point_type is a string that follows the
%                    conventions for POINT_TYPE stated above.
%
% -----------------------
% DATA = COPY(PROB, DATA)

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coco_add_event.m 2839 2015-03-05 17:09:01Z fschild $

%% affected fields in prob
%
%    prob.efunc.events - list of events; events are structures with the
%                        fields
%                        * par  - name of parameter the event is attached
%                                 to
%                        * name - name of event (point type)
%                        * vals - array with values for which to detect an
%                                 event
%                        * type - type of event, one of 'boundary' ('BP'),
%                                 'terminal' ('MX') and 'special point'
%                                 ('SP')

%% check for input argument prob
%  varargin = {[prob], evname|(@evhan, data, [@copy]), [evtype], par, values }

argidx = 1;

if isempty(varargin{argidx}) || isstruct(varargin{argidx})
	prob   = varargin{1};
	argidx = argidx + 1;
else
	prob   = [];
end

if ~isfield(prob, 'efunc') || ~isfield(prob.efunc, 'events')
  prob.efunc.events = [];
end

%% parse event name

event_name    = varargin{argidx};
event_handler = [];
event_handata = [];
event_copy    = [];
argidx        = argidx + 1;

if isa(event_name, 'function_handle')
	event_handler = event_name;
	event_name    = func2str(event_handler);
  event_handata = varargin{argidx};
  argidx        = argidx + 1;
  if isa(varargin{argidx}, 'function_handle')
    event_copy  = varargin{argidx};
    argidx      = argidx + 1;
  end
end

%% parse event type
%  define symbol table for signature constructor functions

evtypes = { ...
	'boundary'      @create_BP_signature ; ...
	'BP'            @create_BP_signature ; ...
	'terminate'     @create_MX_signature ; ...
	'MX'            @create_MX_signature ; ...
	'special point' @create_SP_signature ; ...
	'SP'            @create_SP_signature   ...
};

%  check for evtype argument

idx = find(strcmp(varargin{argidx}, evtypes(:,1)), 1);
if isempty(idx)
	create_signature = @create_SP_signature;
else
	create_signature = evtypes{idx,2};
	argidx           = argidx + 1;
end

%  create signature from remaining arguments
%  assign event name, convert to cell array for later expansion using
%  event.name{[1 1 ... 1]}

event      = create_signature(varargin{argidx:end});
event.name = {event_name};
event.han  = event_handler;
event.data = event_handata;
event.copy = event_copy;

%% append event structure to array prob.events

prob.efunc.events = [ prob.efunc.events; event];


%% parse boundary point signatures
function [ SIG ] = create_BP_signature( varargin )
% EP_SIG = PAR (<|>) VAL | SP_SIG

if nargin==3 && ischar(varargin{2}) && any(varargin{2}=='<>')
	SIG.par  = { varargin{1} };
	SIG.sign =   varargin{2}  ;
	SIG.vals =   varargin{3}  ;

	if numel(SIG.vals) ~= 1
		error('%s: vector of values not allowed in ''<''- or ''>'' signatures', ...
			mfilename);
	end
else
	SIG      = create_SP_signature( varargin{:} );
end

SIG.evlist = 'BP_idx';

%% parse terminal point signatures
function [ SIG ] = create_MX_signature( varargin )
% MX_SIG = PAR [<|>|=] VAL

SIG.evlist = 'MX_idx';

if nargin==3 && ischar(varargin{2}) && any(varargin{2}=='<>=')
	SIG.par  = { varargin{1} };
	SIG.sign =   varargin{2}  ;
	SIG.vals =   varargin{3}  ;
elseif nargin==2
	SIG.par  = { varargin{1} };
	SIG.sign =   '='          ;
	SIG.vals =   varargin{2}  ;
else
	error('%s: wrong number or type of arguments', mfilename);
end

if numel(SIG.vals) ~= 1
	error('%s: number of values must be one for terminal events', mfilename);
end

%% parse special point signatures
function [ SIG ] = create_SP_signature( varargin )
% SP_SIG = PAR [=] VALS | (PAR [=] VAL) ...

if nargin<=0
	error('%s: too few arguments', mfilename);
end

SIG.par    = {};
SIG.vals   = [];
SIG.sign   = '';
SIG.evlist = 'SP_idx';
argidx     = 1;

while argidx<=nargin
	if numel(SIG.par)>=1 && numel(SIG.par)~=numel(SIG.vals)
		error('%s: vector of values not allowed in multi-valued signatures', ...
			mfilename);
	end
	
	SIG.par = [ SIG.par varargin{argidx} ];
	argidx  = argidx + 1;

	if nargin<argidx
		error('%s: too few arguments', mfilename);
	end
	if ischar(varargin{argidx})
		if(varargin{argidx}~='=')
			error('%s: ''%s'': illegal event relation', ...
				mfilename, varargin{argidx});
		end
		argidx = argidx + 1;
	end
	SIG.sign = [SIG.sign '='];

	if nargin<argidx
		error('%s: too few arguments', mfilename);
	end
	SIG.vals = [ SIG.vals(:) ; varargin{argidx}(:) ];
	argidx   = argidx + 1;
end

if numel(SIG.par)>1 && numel(SIG.par)~=numel(SIG.vals)
	error('%s: vector of values not allowed in multi-valued signatures', ...
		mfilename);
end
