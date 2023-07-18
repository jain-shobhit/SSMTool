%========================================================================
% DESCRIPTION:
% Matlab function for continuation of the solution path of a nonlinear
% algebraic system of equations,
%           R(X) = 0,
%           where X = [x;lam],
%           in the interval lam_s <= lam <= lam_e.
% Herein, R and X are real-valued, and R and x have the same dimension.
% The implementation is mainly based on the chapter 'Continuation' in the 
% textbook 'Practical bifurcation and stability analysis' by R. Seydel 
% (2010), and you will find references to equations in that book within
% this code.
%========================================================================
% This file is part of NLvib.
% 
% If you use NLvib, please refer to the book:
%   M. Krack, J. Gross: Harmonic Balance for Nonlinear Vibration
%   Problems. Springer, 2019. https://doi.org/10.1007/978-3-030-14023-6.
% 
% COPYRIGHT AND LICENSING: 
% NLvib Version 1.3 Copyright (C) 2020  Malte Krack  
%										(malte.krack@ila.uni-stuttgart.de) 
%                     					Johann Gross 
%										(johann.gross@ila.uni-stuttgart.de)
%                     					University of Stuttgart
% This program comes with ABSOLUTELY NO WARRANTY. 
% NLvib is free software, you can redistribute and/or modify it under the
% GNU General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% For details on license and warranty, see http://www.gnu.org/licenses
% or gpl-3.0.txt.
%========================================================================
%
% Mandatory variables:
%   x0                  Initial guess
%   fun_residual        Residual function with output R (and dRdX if 
%                       analytical gradients are used) and only input X
%   lam_s,lam_e         Start and end values of parameter
%   ds                  Normalized step size
%
% Optional input variables:
%   varargin{1} = Sopt  Options structure; only specified options are
%                       changed accordingly
%       flag            Fflag whether actual continuation is performed or 
%                       trivial (sequential) continuation is employed 
%                       (default: 1)
%       predictor       Tangent or secant predictors can be specified
%                       ['tangent'|'secant'] (default: 'tangent')
%       parametrization Parametrization constraint
%                       ['arc_length'|'pseudo_arc_length'|'local'| ...
%                           'orthogonal'|'normal'] (default: arc_length)
%       dsmin,dsmax     Minimum, maximum step size (default: ds/5,ds*5)
%       stepadapt       Flag whether step size should be automatically 
%                       adjusted (default: 1; recommended if Sopt.flag = 1)
%       stepmax         Maximum number of steps before termination
%       termination_criterion
%                       Cell array of functions (X) returning logic scalar
%                       1 for termination
%       reversaltolerance
%                       Positive scalar value permitting the path to be
%                       continued below lam_s (if lam_s<lam_e) or beyond
%                       lam_s (if lam_s>lam_e). This could be helpful if
%                       you compute a frequency response e.g. of a system
%                       with softening characteristic and start very close
%                       left to the resonance, since some part of the
%                       solution branch may then be below lam_s. (default:
%                       0.1)
%       jac             Option whether the Jacobian is provided by the
%                       fun_residual; if the Jacobian is only provided
%                       w.r.t. 'x' or not provided at all, it is
%                       approximated using finite differences
%                       ['full'/'on'|'x'|'off'] (default: 'full')
%       Dscale          linear scaling to be applied to vector X
%       dynamicDscale   Flag whether 'Dscale' should be adjusted after
%                       every solution point according to embedded function
%                       'dynamic_Dscale' (default: 0)
%       eps             Tolerance of the norm of the residual in order to
%                       identify failure in convergence
%       noconv_stepcor  In case of failure in convergence, the step size is
%                       either iteratively decreased by factor 2 or
%                       decreased and increased in an alternating manner
%                       ['red'|'redinc']
%       errmax          Maximum number of failure in convergence until the
%                       non-converged solution is accepted or the algorithm
%                       terminates
%       stoponerror     Flag whether algorithm should stop when 'errmax' is
%                       reached (default: 1)
%   varargin{2} = ...
%     fun_postprocess   Cell array of postprocessing function handles e.g.
%                       for stability, sensitivity analyses or amplitude,
%                       energy, response calculation after each successful
%                       iteration
%   varargin{3} = ...
%     Solopt            Options for the nonlinear solver (fsolve)
%
% Output variables:
%   X                   Matrix whose columns are the solution points to
%                       R(X)
%   Solinfo             Structure of solution info
%       NIT             Number of iterations in each step
%       FC              Number of function evaluations in each step
%       IEx             Exit flags of the nonlinear solver in each step
%       ctime           Total computation time
%       FCtotal         Total number of function evaluations
%   Sol                 Structure array whose fields are determined by the
%                       postprocessing functions
function [X,Solinfo,Sol] = ...
    solve_and_continue(x0,fun_residual,lam_s,lam_e,ds,varargin)
%% Handle different input data

% Default continuation options
Sopt = struct('flag',1,'predictor','tangent', ...
    'parametrization','arc_length',...
    'ds',ds,'dsmin',ds/5,'dsmax',ds*5,...
    'stepadapt',1,'stoponerror',1,...
    'reversaltolerance',.1,...
    'stepmax',1e3,'eps',1e-7,'noconv_stepcor','red','errmax',3,...
    'Dscale',ones(length(x0)+1,1),...
    'jac','full','dynamicDscale',0);
if nargin>5 && isstruct(varargin{1})
    % Adapt options if assigned
    tmp = varargin{1};
    new_fieldnames = fieldnames(tmp);
    for ij=1:length(new_fieldnames)
        Sopt.(new_fieldnames{ij}) = ...
            tmp.(new_fieldnames{ij});
    end
end

% Reversal protection
if isfinite(Sopt.reversaltolerance)
    lamlo = min(lam_s,lam_e);
    lamhi = max(lam_s,lam_e);
    lammin = lamlo - abs(Sopt.reversaltolerance*(lamhi-lamlo));
    lammax = lamhi + abs(Sopt.reversaltolerance*(lamhi-lamlo));
    if isfield(Sopt,'termination_criterion')
        Sopt.termination_criterion{end+1} = ...
            @(X) X(end)<lammin || X(end)>lammax;
    else
        Sopt.termination_criterion = {@(X) X(end)<lammin || X(end)>lammax};
    end
end

% Determine continuation direction
if lam_s>lam_e
    dir = -1;
else
    dir = 1;
end

% Adapt variables to scaling
Sopt.ds = Sopt.ds/Sopt.Dscale(end);
Sopt.dsmin = Sopt.dsmin/Sopt.Dscale(end);
Sopt.dsmax = Sopt.dsmax/Sopt.Dscale(end);
lam_s = lam_s/Sopt.Dscale(end);
lam_e = lam_e/Sopt.Dscale(end);
x0 = x0./Sopt.Dscale(1:end-1);

% Postprocessing options
if nargin>6
    fun_postprocess = varargin{2};
else
    fun_postprocess = {};
end

% Solver options
if nargin>7
    Solopt = varargin{3};
else
    Solopt = optimset(optimset(@fsolve),'Display','off',...%'iter',...
        'Jacobian','on','MaxIter',50);%,'DerivativeCheck','on');
end
%% Initialize result vectors
X = zeros(length(x0)+1,Sopt.stepmax);
Sol = struct;
Sol(1:Sopt.stepmax) = struct;
Solinfo.NIT = zeros(Sopt.stepmax,1);
Solinfo.IEx = zeros(Sopt.stepmax,1);
Solinfo.FC = zeros(Sopt.stepmax,1);
Solinfo.ctime = 0;
Solinfo.FCtotal = 0;
%% Initial guess, reference vector, tangent vector and additional row for
% extended Jacobian 'c'
X0 = [x0;lam_s];
zref = [zeros(length(x0),1);1];
z = zref;
Xref = X0;
c = zref;   % c must make rank([J;c']) = n+1, Seydel suggests unit vectors
Xold = X0;
xi = 1;
%% Find initial solution
disp('=================================================================');
disp('NLvib Version 1.3, Copyright (C) 2020 Malte Krack, Johann Gross');
disp('This program comes with ABSOLUTELY NO WARRANTY.');
disp('This is free software, and you are welcome to redistribute');
disp('it under certain conditions, see gpl-3.0.txt.');
disp('=================================================================');
disp('Find initial solution.');
disp('--------------------');

% Relax maximum number of iterations constraint for first solution point
solopt_tmp = Solopt;
Solopt = optimset(solopt_tmp,'MaxIter',500,'Display','iter');

% Solve nonlinear system of equations once for initial solution
if Sopt.flag || isfield(Sopt,'init')&&Sopt.init.flag
    flagtmp = Sopt.flag;
    constrtmp = Sopt.parametrization;
    Sopt.parametrization = 'local';
    if isfield(Sopt,'init')
        yreftmp = Xref;
        zreftmp = zref;
        Xref = X0;
        zref = Sopt.init.zref;
    end
    [XP,~,iEx,output,J] = fsolve(@(X) extended_residual(X, ...
        Xref,zref,fun_residual,Sopt),X0,Solopt);
    Sopt.parametrization = constrtmp;
    X0 = XP;
    if isfield(Sopt,'init')
        Xref = yreftmp;
        zref = zreftmp;
    end
    Sopt.flag = flagtmp;
else
    [XP,~,iEx,output,J] = fsolve(@(X) extended_residual([X;lam_s], ...
        Xref,zref,fun_residual,Sopt),X0(1:end-1),Solopt);
    X0(1:end-1) = XP;
end
if iEx<1
    error('Provided initial guess is not in the basin of attraction.');
end

% Save iteration
X(:,1) = diag(Sopt.Dscale)*X0;
Solinfo.IEx(1) = iEx;
Solinfo.NIT(1) = output.iterations;
Solinfo.FC(1) = output.funcCount;
if ~isempty(fun_postprocess)
    tmp = feval(fun_postprocess,diag(Sopt.Dscale)*X0);
else
    tmp = struct;
end
new_fieldnames = fieldnames(tmp);
for ij=1:length(new_fieldnames)
    Sol(1).(new_fieldnames{ij}) = ...
        tmp.(new_fieldnames{ij});
end
disp(['Continuation at ' num2str(X0(end)*Sopt.Dscale(end),'%.4f') ...
    ', step size ' num2str(Sopt.ds*Sopt.Dscale(end)) '.']);

% Reset maximum number of iterations
Solopt = solopt_tmp;
%% Start continuation
disp('--------------------');
disp('Start continuation.');
disp('--------------------');
tic
ierr = 0;
istep = 2;
while istep<=Sopt.stepmax
    %% RESCALE (dynamically adjust scaling, if requested)
    if Sopt.dynamicDscale
        if istep == 2
            % Store initial Dscale for dynamic adoption
            Sopt.Dscale0 = Sopt.Dscale;
        end
        %% Adjust variable scaling values (diagonal elements of diagonal
        % scaling matrix) so that scaled variables have value ~1. But avoid
        % weighing very small values too much by setting as minimum dynamic
        % scaling the initial one (Dscale0).
        Dscaleold = Sopt.Dscale;
        Sopt.Dscale(1:end-1) = max(abs(X0(1:end-1).*Dscaleold(1:end-1)),...
            Sopt.Dscale0(1:end-1));
        %% Update scaled variable values, Jacobian, reference tangent
        Xref    = Xref.*(Dscaleold./Sopt.Dscale);
        X0      = X0.*(Dscaleold./Sopt.Dscale);
        Xold    = Xold.*(Dscaleold./Sopt.Dscale);
        if Sopt.flag
            J       = J*diag(1./Dscaleold);
            zref    = zref.*Dscaleold;
        end
    end
    %% PREDICT
    
    % Determine predictor direction (in the unscaled system)
    if Sopt.flag
        switch Sopt.predictor
            case 'tangent'
                [~,kk] = sort(abs(zref./max(abs(Xref),1e-4)),...
                    1,'descend');
                % Temporarily switch off warning
                warning('off','MATLAB:nearlySingularMatrix');
                warning('off','MATLAB:singularMatrix');
                for ik=1:length(kk)
                    k = kk(ik);
                    c = zeros(length(X0),1); c(k) = 1;
                    % Determine unit tangent to the solution
                    % path (Eq. 4.8)
                    ztmp = [J(1:end-1,:);c']\...
                        [zeros(size(J,1)-1,1);1];
                    if ~any(isnan(ztmp))
                        % Successful!
                        break;
                    end
                end
                % Switch warning on again
                warning('on','MATLAB:nearlySingularMatrix');
                warning('on','MATLAB:singularMatrix');
            case 'secant'
                if istep>2
                    % Secant predictor
                    ztmp = X0-Xold;
                else
                    % Take tangent step at first iteration
                    ztmp = [J(1:end-1,:);c']\...
                        [zeros(size(J,1)-1,1);1];
                end
            otherwise
                error('Unknown predictor specification.');
        end
        if any(isnan(ztmp))
            error('Could not determine predictor direction.');
        end
        % Apply linear scaling to tangent
        ztmp = ztmp./Sopt.Dscale;
        z = ztmp/norm(ztmp);
    else
        % Trivial tangent, already initialized
    end
    
    % Take step
    XP = X0+dir*Sopt.ds*z;
    
    % Ensure forward stepping along solution path
    if Sopt.flag
        if (istep > 2) && ...
                (transpose(X0-Xold)*(dir*Sopt.ds*z) < 0);
            XP = X0-dir*Sopt.ds*z;
        end
    end
    %% CORRECT
    if Sopt.flag
        % Determine reference data for parametrization equation
        switch Sopt.parametrization
            case 'local'
                % Is supposed to work well in conjunction with secant
                % predictor
                zref = zeros(size(X0));
                Xref = zeros(size(X0));
                if istep>2
                    % Local parametrization Eq. (4.18-19)
                    
                    % Find variables with maximum relative change at curent
                    % point
                    [~,kind] = sort(abs(XP-X0)./max(abs(X0),1e-4));
                    
                    % In case of error, take different index for local
                    % parametrization
                    k = kind(end-min(ierr,length(kind)-1));
                    Xref(k) = XP(k);
                else
                    k = length(X0);
                    Xref(k) = XP(k);
                end
                zref(k) = 1;
            case {'orthogonal','normal'}
                % For orthogonal or normal parametrization, the predicted 
                % solution is relevant
                Xref = XP;
                zref = z;
            case {'arc_length','pseudo_arc_length'}
                % For arc-length-like parametrization, the previous
                % solution is relevant
                Xref = X0;
                zref = z;          
            otherwise
                error(['''',Sopt.parametrization,'''', ...
                    ' is an invalid specification'...
                    ' of path continuation' ...
                    ' parametrization constraint.']);
        end
        % Solve extended nonlinear system of equations
        [Xtmp,Rext,iEx,output,Jtmp] = ...
            fsolve(@(X) extended_residual( ...
            X,Xref,zref,fun_residual,Sopt),XP,Solopt);
    else
        % Solve extended nonlinear system of equations
        [Xtmp,Rext,iEx,output,Jtmp] = ...
            fsolve(@(X) extended_residual( ...
            [X;XP(end)],Xref,zref,fun_residual,Sopt),XP(1:end-1),Solopt);
    end
    %% ACCEPT/REJECT POINT
    if iEx<1 && sqrt(Rext'*Rext)>Sopt.eps
		% Reject point
        if  (Sopt.ds == Sopt.dsmin) || (ierr>=Sopt.errmax)
            if Sopt.stoponerror
                % Delete unused fields and break
                X(:,istep:end) = [];
                Solinfo.NIT(istep:end) = [];
                Solinfo.IEx(istep:end) = [];
                Solinfo.FC(istep:end) = [];
                disp(['No convergence, stopping ' ...
                    'continuation.']);
                break;
            else
                disp(['No convergence, proceed computation with' ...
                    ' unconverged solution norm(Rres)=' num2str(norm(Rext))]);
                ierr = 0;
            end
        else
            ierr = ierr+1;
            if ierr==1
                ds_ref = Sopt.ds;
                xiref = xi;
            end
            % Adapt stepsize
            switch Sopt.noconv_stepcor
                case 'red'
                    fact = ierr+1;
                    Sopt.ds = max(ds_ref/2^fact,Sopt.dsmin);
                    disp(['No convergence, reduce stepsize to ds = ' ...
                        num2str(Sopt.ds*Sopt.Dscale(end))]);
                case 'redinc'
                    fact = fix((ierr+1)/2)+1;
                    if mod(ierr,2)
                        Sopt.ds = max(ds_ref/2^fact,Sopt.dsmin);
                        disp(['No convergence, reduce stepsize' ...
                            ' to ds = ' num2str(Sopt.ds*Sopt.Dscale(end))]);
                    else
                        Sopt.ds = min(ds_ref*2^fact,Sopt.dsmax);
                        disp(['No convergence, increase stepsize' ...
                            ' to ds = ' num2str(Sopt.ds*Sopt.Dscale(end))]);
                    end
                case 'none'
                otherwise
                    error(['Unknown step correction specifier' ...
                        ' in case of fail in convergence.']);
            end
            xi = xiref*Sopt.ds/ds_ref;
            if ierr<Sopt.errmax
                continue;
            else
                ierr = 0;
            end
        end
    else
        ierr = 0;
    end
    %% POSTPROCESS SOLUTION POINT
    
    % Save previous solution
    Xold = X0;
    % Update solution vector
    if Sopt.flag
        X0 = Xtmp;
        J = Jtmp;
    else
        X0 = [Xtmp;XP(end)];
    end
    
    % Save solution point
    X(:,istep) = diag(Sopt.Dscale)*X0;
    Solinfo.IEx(istep) = iEx;
    Solinfo.NIT(istep) = output.iterations;
    Solinfo.FC(istep)  = output.funcCount;
    
    % Call postprocessing function
    if ~isempty(fun_postprocess)
        Sol(istep) = feval(fun_postprocess,diag(Sopt.Dscale)*X0);
    end
    %% TERMINATE/PROCEED CONTINUATION
    
    % Evaluate termination criteria
    if isfield(Sopt,'termination_criterion')
        term = zeros(length(Sopt.termination_criterion),1);
        for iterm=1:length(Sopt.termination_criterion)
            term(iterm) = feval(Sopt.termination_criterion{iterm},...
                X(:,istep));
        end
    else
        term = 0;
    end
    
    % Look if end of interval is reached or any of the termination criteria
    % is met
    if ( ((X0(end) >= lam_e)&&(dir==1)) || ...
            ((X0(end) <= lam_e)&&(dir==-1)) ) || any(term)
        % Delete unused fields and break
        X(:,istep+1:end) = [];
        Sol(istep+1:end) = [];
        Solinfo.NIT(istep+1:end) = [];
        Solinfo.IEx(istep+1:end) = [];
        Solinfo.FC(istep+1:end) = [];
        if any(term)
            disp(['Terminating continuation since at least one of the ' ...
                'termination criteria is met.']);
        else
            disp(['Terminating continuation since parameter end value ' ...
                'is reached.']);
        end
        break;
    else
        disp(['Continuation at ' num2str(X0(end)*Sopt.Dscale(end),'%.4f') ...
            ', step size ' num2str(Sopt.ds*Sopt.Dscale(end)) '.']);
    end
    
    if istep==Sopt.stepmax
        disp(['Terminating continuation since maximum number of ' ...
            'solution points ' num2str(Sopt.stepmax) ' is reached.']);
    end
    %% ADAPT STEP LENGTH
    if istep>2 && Sopt.stepadapt
        % Increase/Decrease step size if number of solver
        % iterations is smaller/larger than specified number
        if all(Solinfo.NIT(istep-1:istep)<5) || ...
                (istep>4&&all(Solinfo.NIT(istep-4:istep)<6))
            xi = 2;
        elseif all(Solinfo.NIT(istep-1:istep)>9)
            xi = 0.5;
        else
            xi = 1;
        end
        dstmp = xi*Sopt.ds;
        
        % Impose bound constraints on step size
        Sopt.ds = max(dstmp,Sopt.dsmin);
        Sopt.ds = min(Sopt.ds,Sopt.dsmax);
    end
    
    % Increment loop count
    istep = istep+1;
end
%% Output computational effort
Solinfo.ctime = toc;
Solinfo.FCtotal = sum(Solinfo.FC);
disp('--------------------');
disp('COMPUTATIONAL EFFORT:');
disp(['Total elapsed time (toc) is ',num2str(Solinfo.ctime,'%16.1f'),' s']);
disp(['Total number of function evaluations is ', ...
    num2str(Solinfo.FCtotal)]);
end
function [Rext,dRextdX] = extended_residual(X,Xref,zref,...
    fun_residual,Sopt)
%% Evaluation of the residual function and its derivative
switch Sopt.jac
    case {'full','on'}
        [R,dRdX] = feval(fun_residual,diag(Sopt.Dscale)*X);
    case 'x'
        [R,dRdx] = feval(fun_residual,diag(Sopt.Dscale)*X);
        
        % Approximate dfdlam using finite differences
        dlam = max(Sopt.eps*abs(X(end)),Sopt.eps);
        Xtmp = X; Xtmp(end) = X(end)+dlam;
        dlam = Xtmp(end)-X(end);
        Rp = feval(fun_residual,diag(Sopt.Dscale)*Xtmp);
        dRdlam = (Rp-R)/dlam/Sopt.Dscale(end);
        dRdX = [dRdx dRdlam];
    otherwise
        R = feval(fun_residual,diag(Sopt.Dscale)*X);
        Smyopt = struct('epsrel',Sopt.eps,'epsabs',Sopt.eps, ...
            'ikeydx',1,'ikeyfd',1);
        dRdX = finite_difference_jacobian(fun_residual,...
            diag(Sopt.Dscale)*X,Smyopt);
end
%% Evaluation of the parametrization constraint equation and its derivative
if Sopt.flag
    switch Sopt.parametrization
        case {'arc_length'}
            % iteration on a normal plane, perpendicular to tangent
            p = transpose((X-Xref))*(X-Xref)-Sopt.ds^2;
            dpdX = 2*transpose(X-Xref);
        case 'pseudo_arc_length'
            p = Sopt.pseudoxi*...
                transpose((X(1:end-1)-...
                Xref(1:end-1)))*(X(1:end-1)-...
                Xref(1:end-1)) + (1-Sopt.pseudoxi)*...
                (X(end)-Xref(end))^2 - ...
                Sopt.ds^2;
            dpdX = [Sopt.pseudoxi*...
                2*transpose(X(1:end-1)-Xref(1:end-1)) ...
                2*(1-Sopt.pseudoxi)*...
                (X(end)-Xref(end))] ;
        case {'local','orthogonal'}
            % Solution on hyperplane through Xref, normal to zref
            % Eq. (4.22)
            p       = transpose(zref)*(X-Xref);
            dpdX    = transpose(zref);
        case {'normal'}
            [~,kk] = sort(abs(zref./max(abs(X),1e-4)),...
                1,'descend');
            % Temporarily switch off warning
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:singularMatrix');
            for ik=1:length(kk)
                k = kk(ik);
                c = zeros(length(X),1); c(k) = 1;
                % Determine unit tangent to the solution
                % path (Eq. 4.8)
                ztmp = [dRdX;c']\...
                    [zeros(size(dRdX,1),1);1];
                if ~any(isnan(ztmp))
                    % Erfolgreich!
                    break;
                end
            end
            % Switch warning on again
            warning('on','MATLAB:nearlySingularMatrix');
            warning('on','MATLAB:singularMatrix');
            z       = ztmp/norm(ztmp);
            p       = z'*(X-Xref);
            dpdX    = transpose(z);
        otherwise
            error('Invalid specification of path continuation constraint');
    end
end
%% Assembly of extended residual and derivative
if Sopt.flag
    Rext = [R;p];
    dRextdX = [dRdX*diag(Sopt.Dscale);dpdX];
else
    Rext = R;
    dRextdX = dRdX(:,1:end-1)*diag(Sopt.Dscale(1:end-1));
end
end
function J = finite_difference_jacobian(Hfuncname, x0, Smyopt, varargin)
%========================================================================
% DESCRIPTION:
% Function for the calculation of a Jacobian matrix. The derivatives of the
% Jacobian matrix are approximated by forward finite differences, by back-
% ward finite differences or by central finite differences. Different ways
% are possible for the calculation of the step size needed for the finite
% differences approximation of the derivatives. The algorithm is based on
% the FORTRAN77 subroutine 'fdjac' presented in the book 'Numerical Recipes
% in FORTRAN77' (chapter 9, section 7, pages 376-386), Cambridge University
% Press, 1992.
%========================================================================
%
% Called functions:
%   Hfuncname (function handle)
%
% Input data dictionary:
%   Hfuncname           FunctionHandle  Function handle for the function
%   Smyopt.epsrel       Struct.Rscalar  Maximum relative error
%   Smyopt.epsabs       Struct.Rscalar  Machine precision tolerance
%   Smyopt.ikeydx       Struct.Iscalar  Key for the step size calculation
%   Smyopt.ikeyfd       Struct.Iscalar  Key for the finite differences
%                                       method
%   x0                  Rarray          Input variable vector
%   varargin            Cell Array      Optional parameters for the
%                                       function
%
% Output data dictionary:
%   J                   Rarray          Computed Jacobian matrix
%% Start of the calculation of the Jacobian matrix based on finite
% differences
for ii = 1:length(x0)
    
    %% Calculation of the actual step size for the finite difference method
    switch (Smyopt.ikeydx)
        
        % First method for the calculation of the finite difference step
        % size
        case (1)
            dx = Smyopt.epsrel*abs(x0(ii));
            if (dx == 0)
                dx = Smyopt.epsrel;
            end
            
            % Second method for the calculation of the finite difference
            % step size
        case (2)
            dx = sqrt(Smyopt.epsabs)*(1+abs(x0(ii)));
            
            % Third method for the calculation of the finite difference
            % step size
        case (3)
            dx = Smyopt.epsrel*abs(x0(ii));
            dx = max(dx,Smyopt.epsabs);
            
    end
    %% Jacobian matrix computation based on different finite difference methods
    switch (Smyopt.ikeyfd)
        
        % Calculation based on a forward finite difference approximation
        case (1)
            
            % Initialize the Jacobian matrix and the function value
            if (ii == 1)
                f0 = feval(Hfuncname,x0,varargin{:});
                J = zeros(length(f0),length(x0));
            end
            
            % Udpate of the input variable vector
            xp = x0;
            xp(ii) = x0(ii)+dx;
            
            % Recompute the finite difference step size to reduce the
            % finite precision error
            dx = xp(ii)-x0(ii);
            
            % Calculation of the actual column of the Jacobian matrix
            fp = feval(Hfuncname,xp,varargin{:});
            J(:,ii) = (fp-f0)/dx;
            %% Calculation based on a backward finite difference approximation
        case (2)
            
            % Initialize the Jacobian matrix and the function value
            if (ii == 1)
                f0 = feval(Hfuncname,x0,varargin{:});
                J = zeros(length(f0),length(x0));
            end
            
            % Udpate of the input variable vector
            xm = x0;
            xm(ii) = x0(ii)-dx;
            
            % Recompute the finite difference step size to reduce the
            % finite precision error
            dx = x0(ii)-xm(ii);
            
            % Evaluate the function value at the updated input variable
            % vector
            fm = feval(Hfuncname,xm,varargin{:});
            
            % Calculation of the actual column of the Jacobian matrix
            J(:,ii) = (f0-fm)/dx;
            %% Calculation based on a central finite difference approximation
        case (3)
            
            % Udpate of the input variable vector
            xp = x0;
            xp(ii) = x0(ii)+dx;
            xm = x0;
            xm(ii) = x0(ii)-dx;
            
            % Recompute the finite difference step size to reduce the
            % finite precision error
            dx = (xp(ii)-xm(ii))/2;
            
            % Evaluate the function value at the updated input variable
            % vector
            fp = feval(Hfuncname,xp,varargin{:});
            fm = feval(Hfuncname,xm,varargin{:});
            
            % Initialize the Jacobian matrix
            if (ii == 1)
                J = zeros(length(fp),length(x0));
            end
            
            % Calculation of the actual column of the Jacobian matrix
            J(:,ii) = (fp-fm)/2/dx;
    end
end
end