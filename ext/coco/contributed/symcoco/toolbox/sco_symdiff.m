function [df,dx]=sco_symdiff(f,x,varargin)
%% directional derviatives of fs wrt xxs and ps
%
%% Inputs:
%
% * f: m x 1 array of symbolic expressions using xxs as states and ps as
% parameters
% * x: nx1 array of inputs
%
%% Optional inputs (name-value pairs) 
%
% * maxorder: 2 
% * append='_dev': names used for deviation directions, make sure these do
% not clash with original names in x.
% * deviation_name: default |h_devsmall| deviation variable (will be set to
% 0 after differentiation)
%
%% Outputs:
%
% * df: cell array of length maxorder where |df{k}| contains the the kth
% directional derivative.
% $Id$
%%
default={'maxorder',2,'dev_append','_dev','deviation_name','h_devsmall'};
options=sco_set_options(default,varargin,'pass_on');
n=length(x);
%% define deviations
xnames=sco_names_from_sym(x);
xdevnames=cell_name_append(xnames,options.dev_append);
xnames=xnames(:).';
xdevnames=xdevnames(:).';
hname=options.deviation_name;
if ~isempty(intersect(xnames,xdevnames)) || ismember(hname,[xnames,xdevnames])
    error('symdiff:names',...
        'symdiff: name clash: choose different option for append or deviation_name\n%s\n%s\n%s\n',...
        evalc('celldisp(xnames,''variables'')'),...
        evalc('celldisp(xdevnames,''deviations'')'),...
        evalc('celldisp(xnames,''hname'')'));
end
dx=reshape(cell2sym(xdevnames),[n,1]);
hdev=sym(options.deviation_name);
%% add deviations to xx and par
fdev=subs(f,x,x+hdev*dx);
%% differentiate i times
df=cell(1,options.maxorder);
for i=1:options.maxorder
    hi=repmat({hdev},1,i);
    deriv=diff(fdev,hi{:});
    df{i}=subs(deriv,hdev,0);
end
end
%% append letter to cell array of names
function xnap=cell_name_append(xn,app)
xnap=cell(size(xn));
for i=1:numel(xn)
    xnap{i}=[xn{i},app];
end
end

