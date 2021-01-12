function str=sco_symcode(df,x,dx,varargin)
%% convert symbolic expressions to matlab code using matlabFunction (uses temporary files!)
%
% output is string containing function body with name chooseable by
% optional input funcname.
%
% set optional input scalar to true to return every component as separate
%
% $Id: dde_symcode.m 193 2017-06-12 10:31:06Z jansieber $
%%
%#ok<*AGROW>
default={'funcname','sys','keeptemp',false,'outname','out'};
options=sco_set_options(default,varargin,'pass_on');
%% generate code with matlabFunction
% matlabFunction can only write to files, so generate temporary files where
% we save the function then read the file as string, concatenated 
nl=sprintf('\n'); %#ok<SPRINTFN>
str='';
if ~sco_isoctave()
    folder=tempname;
    gotfolder=mkdir(folder);
    if ~gotfolder
        error('sco_symcode:perm','sco_symcode: could not create temp folder %s',folder);
    end
end
%% break up output array into sequence of outputs 
% to enable scalar expansion after function call
vars=[x(:).',dx(:).'];
vnames=sco_names_from_sym(vars);
nf=length(df{1});
if sco_isoctave()
    %% octave sym cannot handle setting output names
    options.out='out';
end
outnames=arrayfun(@(i)sprintf('%s_%d',options.outname,i),1:nf,'uniformoutput',false);
if sco_isoctave()
    outargs={'outputs',outnames};
else
    outargs={};
end
%%
if ~isempty(intersect(vnames,outnames))
    error('sco_symcode:names',...
        'sco_symcode: name clash between output names and variables');
end
%% generate code in files tempname/funcname_(order-1).m and re-read into string str
for ic=1:length(df)
    dfcell=num2cell(df{ic});
    fname=sprintf('%s_%d',options.funcname,ic-1);
    if sco_isoctave()
        strnew=sco_octave_code(dfcell,fname,vars);
    else
        % write code to temporary file
        filename=fullfile(folder,[fname,'.m']);
        w_orig=warning;
        warning('off','symbolic:generate:FunctionNotVerifiedToBeValid');
        matlabFunction(dfcell{:},'file',filename,'vars',vars,outargs{:});
        warning(w_orig);
        % read code back into string str
        fid=fopen(filename,'r');
        strnew=fread(fid,inf);
        fclose(fid);
    end
    str=[str,nl,char(strnew(:)'),nl];
end
%% remove folder (unless optionally prevented)
if ~sco_isoctave && ~options.keeptemp
    rmdir(folder,'s')
end
end
%% Simplified version of function_handle from OctSymPy
%
% returns only code string
% adapted from original functiona_handle (Copyright (C) 2014-2016 Colin
% B. Macdonald)
function code = sco_octave_code(f,fcnname,inputs)
%% command as in function_handle, however with pre-set arguments
cmd = { '(expr,fcnname,in_vars) = _ins' ...
    'from sympy.utilities.codegen import codegen' ...
    'try:' ...
    ['    out = codegen((fcnname,expr), "' 'octave' ...
    '", fcnname, header=1' ...
    ', argument_sequence=in_vars)'] ...
    'except ValueError as e:' ...
    '    return (False, str(e))' ...
    'return (True, out)' };
[worked, out] = python_cmd (cmd, f, fcnname, inputs);
if (~worked)
    if (strcmp(out, 'Language ''octave'' is not supported.'))
        error('function_handle: your SymPy has no octave codegen');
    else
        out
        error('function_handle: Some other error from SymPy code gen?  file a bug!');
    end
end
code = out{1}{2};
end
