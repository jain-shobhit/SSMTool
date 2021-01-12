function [funcstr,derivatives]=sco_sym2funcs(f,args,names,varargin)
%% create matlab code for r.h.s (and delays if needed) from symbolic expressions
%
%   function [funcstr,derivatives]=sco_sym2funcs(f,args,names,...) 
%
% It converts the input f (an array of symbolic expressions) into a
% right-hand side function and its derivatives that can be used for
% the coll toolbox. 
%
%% Inputs:
%
% * |f|:  n x 1 array of symbolic expressions, defining the right-hand side
% * |args|: cell array of length(args{i})x1 long input arguments
% * |names|: list of names used to extract derivatives
%
%% Common optional name-value input pairs
% * |'vector'|: logical array of same size as args, a hint for postprocessing,
% whether this argument should be treated as scalar or vector
% * |'write'| (default |true|): write output to file?
% * |'filename'| (default |'sys'|): results are written to function file
% with this name.
% * |'maxorder'| (default=max 2): maximal order of derivatives to be computed.
%
%% Outputs (often not needed)
%
% * |fstr|: strong containing the code for writing the function.
% * |derivatives|: symbolic expressions for all dervatives (a structure,
% containing the fields |df|, |xx|, |parameter|, |dx| and |dp|. |df|
% contains the expressions, |xx|, |parameter|, |dx| and |dp| the symbols
% used in these expressions. 
%
% The routine will create a vectorized matlab
% function |y=f(args{:})| from the expressions in |fs|, and its directional
% derivatives up to a maximum order
%% Warning
% The file will write to temporary files, since matlabFunction can only
% write to files.
%
% $Id$
%%
default={'vector',true(1,length(args)),'write',true,...
    'filename','sys','extension','rhs'};
[options,pass_on]=sco_set_options(default,varargin,'pass_on');
%% convert possibly cell arrays safely into sym column vectors
f=sco_sym_from_cell(f);
nargs=length(args);
len=NaN(1,nargs);
for i=nargs:-1:1
    args{i}=sco_sym_from_cell(args{i});
    args{i}=args{i}(:);
    len(i)=length(args{i});
end
f=f(:);
nf=length(f);
rg_bounds=[0,cumsum(len)];
rg_bounds=[rg_bounds(1:end-1)+1;rg_bounds(2:end)]';
rgstr=create_str(names,rg_bounds);
szstr=create_str(names,len(:));
vecstr=create_str(names,options.vector(:));
x=vertcat(args{:});
%% extract function name from (optional) filename
[pth,funcname,ext]=fileparts(options.filename);
%% differentiate
[df,dx]=sco_symdiff(f,x,pass_on{:});
df=[{f},df];
%% generate code
str=sco_symcode(df,x,dx,'funcname',[funcname,'_',options.extension],pass_on{:});
%[fstr{1},df{1},dx{1}]=sco_symdericode(f,x,[funcname,'_',options.extension],pass_on{:});
maxorder=length(df)-1;
%str=[fstr{:}];
derivatives=struct('df',df,'x',x,'dx',dx);
%% octave's output ends with 'end', matlab's does not
if sco_isoctave()
    function_end='end';
else
    function_end='';
end
%% create full function (still string)
nl=sprintf('\n'); %#ok<SPRINTFN>
header=sprintf('function varargout=%s(action,varargin)',funcname);
comment=[...
    '%% Automatically generated with matlabFunction',nl,...
    '%#ok<*DEFNU,*INUSD,*INUSL>',nl];
body=[...
    'switch action',nl,...
    '  case ''nargs''',nl,...
    '   varargout{1}=',num2str(nargs),';',nl,...
    '   return',nl,... 
    '  case ''nout''',nl,...
    '   varargout{1}=',num2str(nf),';',nl,...
    '   return',nl,... 
    '  case ''argrange''',nl,...
    '   varargout{1}=',rgstr,nl,...
    '   return',nl,... 
    '  case ''argsize''',nl,...
    '   varargout{1}=',szstr,nl,...
    '   return',nl,... 
    '  case ''vector''',nl,...
    '   varargout{1}=',vecstr,nl,...
    '   return',nl,... 
    '  case ''extension''',nl,...
    '   varargout{1}=''',options.extension,''';',nl,...
    '   return',nl,... 
    '  case ''maxorder''',nl,...
    '   varargout{1}=',num2str(maxorder),';',nl,...
    '   return',nl,... 
    'end',nl,...
    'nout=',num2str(nf),';',nl,...
    'order=varargin{1};',nl,...
    'f=str2func(sprintf(''',funcname,'_%s_%d'',action,order));',nl,...
    'varargout=cell(nout,1);',nl,...
    '[varargout{:}]=f(varargin{2:end});',nl,...
    function_end,nl,...
    nl];
funcstr=[header,nl,comment,nl,body,nl,str];
%% write to file if requested
if options.write
    if isempty(ext)
        ext='.m';
    end
    filename=fullfile(pth,[funcname,ext]);
    fid=fopen(filename,'w');
    if fid<=0
        error('sco_sym2funcs:perm','sco_sym2funcs: could not create function file %s',filename);
    end
    fwrite(fid,funcstr);
    fclose(fid);
end
end
%% create structure containing function arguments
function rgstr=create_str(names,rgval)
rgstr='struct(';
for i=1:length(names)
    rgstr=[rgstr,'''',names{i},''',']; %#ok<AGROW>
    if size(rgval,2)==2
        rgstr=[rgstr,num2str(rgval(i,1)),':',num2str(rgval(i,2)),',']; %#ok<AGROW>
    else
        rgstr=[rgstr,num2str(rgval(i,1)),',']; %#ok<AGROW>
    end
end
rgstr=[rgstr(1:end-1),');'];
end
%% Convert possibly cells to sym arrays safely
% This wrapper converts if necessary, but does nothing otherwise.
function xa=sco_sym_from_cell(xc)
if iscell(xc)
    xa=cell2sym(xc);
elseif isa(xc,'sym')
    xa=xc;
else
    error('sco_sym_from_cell:input',[...
        'sco_sym_from_cell: input is neither cell (of strings) nor sym array',...
        'but of type %s'],class(xc));
end
end
