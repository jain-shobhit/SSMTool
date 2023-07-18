%========================================================================
% DESCRIPTION:
% Matlab function applying the Newmark average constant acceleration time 
% step integration method (implicit, unconditionally stable, 
% low/no numerical damping for simple cases, cf. textbooks) to the problem
% 
%       system.M * \ddot q + system.D * \dot q + system.K * q + ...
%                   f_nl(q,\dot q) - f_ex(t) = 0            (1)
% 
%          with initial values      y(0) = [q(0);\dot q(0)] = ys
%          within time interval     0 < t <= Np*2*pi/Om.
% 
% The method uses 'Ntd' equidistant time steps per period.
% For numerical reasons, the problem is time-normalized with respect to the
% oscillation period (tau = Om*t, dq/dt = dq/dtau*Om, and so on).
% 
% The function sets up the vector 'R=ye-ys', i.e. the deviation between 
% initial and end values, which is the residual vector of the shooting
% method. Of course, the function can be also used as integrator.
% 
% The relations between displacement, velocities and accelerations at start
% 'S' and end 'E' of a time step for the above mentioned algorithm are
% 
%       \dot qE = \dot qS + (\ddot qS + \ddot qE)/2 * dt,   (2a)
%            qE = qS + (\dot qS + \dot qE)/2 * dt.          (2b)
% 
% When we evaluate Eq.(1) at tE and substitute Eq.(2a-b), we can derive the
% implicit nonlinear (if f_nl is nonlinear) algebraic equation for each
% time step,
% 
%           S*qE + f_nl(qE) - b(tE,qE) = 0                  (3)
% 
%       with          S = 4/dt^2*M + 2/dt*D + K 
%       and    b(tE,qE) = f_ex(tE) + ...
%                           M*(4/dt^2*qS + 4/dt*\dot qS + \ddot qS) + ...
%                           D*(2/dt*qS + \dot qS),
% which we solve for qE with a fixed-point iteration scheme.
% 
% Let us write \dot q as u in this script!
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
function [R,dR,ye,Y,dye_dys] = shooting_residual(X,system,...
    Ntd,Np,analysis,qscl,fscl,inorm)
%% Handle analysis type
if nargin<6 || isempty(qscl)
    qscl = 1;
end
if nargin<7 || isempty(fscl)
    fscl = mean(diag(system.K))*qscl;
end
n = system.n;
switch lower(analysis)
    case {'frf','frequency response'}
        % Frequency response analysis: X = [ys;Om]
        
        % Initial displacement and velocity
        ys = X(1:2*n);
        q = ys(1:n); u = ys(n+1:2*n);
        dX = eye(length(X));
        dys = dX(1:2*n,:);
        dq = dX(1:n,:); du = dX(n+1:2*n,:);
        
        % Excitation 'excitation' is the fundamental harmonic of the external
        % forcing
        fex = @(tau) real(system.Fex1*exp(1i*tau));
        
        % Excitation frequency
        Om = X(end);
        dOm = zeros(1,length(X)); dOm(end) = 1;
        
        % Mass, damping and stiffness matrices of the time-normalized problem
        M = Om^2*system.M;  dM_dOm = 2*Om*system.M;
        D = Om*system.D;    dD_dOm = system.D;          dD_ddel = 0*dD_dOm;
        K = system.K;
        
        % Damping 'del' does not exist in this case, set its derivative to
        % zero
        ddel = zeros(1,size(X,1));
    case {'nma','nonlinear modal analysis'}
        % Nonlinear modal analysis:  X = [qs_/a;us_/(a*om);om;D;log10a]
        dX = eye(length(X));
        
        % Modal amplitude
        a = exp(log(10)*X(end));
        da = log(10)*exp(log(10)*X(end))*dX(end,:);
        
        % Oscillation frequency
        Om = X(end-2); dOm = dX(end-2,:);
        
        % Modal damping ratio
        del = X(end-1); ddel = dX(end-1,:);
        
        % Mass, damping and stiffness matrices of the time-normalized
        % problem
        M = Om^2*system.M;
        dM_dOm = 2*Om*system.M;
        D = Om*system.D - 2*del*Om^2*system.M;
        dD_dOm = system.D - 4*del*Om*system.M;
        dD_ddel = -2*Om^2*system.M;
        K = system.K;
        
        % Initial displacement and velocity
        nnorm = setdiff(1:2*n,[inorm inorm+n]);
        ys = zeros(2*n,1);          dys = zeros(size(ys,1),size(X,1));
        ys(inorm) = a;              dys(inorm,:) = da;
        ys(nnorm) = a*X(1:end-3);   dys(nnorm,:) = a*dX(1:end-3,:) + X(1:end-3)*da;        
        q = ys(1:n);                dq = dys(1:n,:);
        u = ys(n+1:2*n);            du = dys(n+1:2*n,:);
        
        % No forcing
        fex = @(tau) 0;
    otherwise
        error(['Unknown analysis type ' analysis '.']);
end
%% Setup parameters of the integration

% Time step size (accordingly normalized)
dt = 2*pi/Ntd;

% Matrix involved in the implicit one-step problem (constant if we do not
% change the time step size or the algorithm)
S = 4/dt^2*M+2/dt*D+K;
dS_dOm = 4/dt^2*dM_dOm+2/dt*dD_dOm;
dS_ddel = 2/dt*dD_ddel;

% Initialize output, if requested
if nargout>3
    Y = zeros(Ntd,2*n);
end
%% Calculate initial accelerations (account for time-normalization: Om*u is
% the non-normalized velocity)
[fnl,dfnl_dq,dfnl_du] = nonlinear_forces(q,Om*u,system.nonlinear_elements);
acc = M\(fex(0)-D*u-K*q-fnl);
dacc = -M\( D*du+dD_dOm*u*dOm+dD_ddel*u*ddel + K*dq + ...
    dfnl_dq*dq+dfnl_du*(Om*du+u*dOm)) - ...
    2*acc/Om*dOm;
%% Loop over time steps
for j=1:Np*Ntd
    %% Calculate displacement at next time instant
    
    % Determine predicted displacement
    qE_pred = q + u*dt;
    
    % Evaluate constant vector 'b' involved in the implicit one-step
    % problem
    b = fex(j*dt) + M*( 4/dt^2*q + 4/dt*u + acc ) + D*( 2/dt*q + u );
    db = M*(4/dt^2*dq+4/dt*du+dacc) + D*(2/dt*dq+du) + ...
        dM_dOm*( 4/dt^2*q + 4/dt*u + acc )*dOm + ...
        dD_dOm*( 2/dt*q + u )*dOm + dD_ddel*( 2/dt*q + u )*ddel;

    % Solve for nonlinear displacement at end of time step
    [qE,dqE] = ...
        solve_qE(qE_pred,S,b,system.nonlinear_elements,q,u,Om,dt,...
        dS_dOm,dS_ddel,db,dq,du,dOm,ddel,fscl);
    
    % Evaluate velocity and acceleration at next time level
    uE = 2/dt*(qE-q)-u;
    accE = 4/dt^2*(qE-q)-4/dt*u-acc;
    duE = 2/dt*(dqE-dq)-du;
    daccE = 4/dt^2*(dqE-dq)-4/dt*du-dacc;
    
    % Save displacement and velocity of current time level, if requested
    if nargout>3
        Y(j,:) = [q' u'];
    end
    
    % Update current values of displacement, velocity and acceleration
    q = qE; u = uE; acc = accE;
    dq = dqE; du = duE; dacc = daccE;
end
%% Calculate residual and Jacobian
ye = [q;u];
R = (ye-ys)/qscl;
dye_dys = [dq(:,1:2*n);du(:,1:2*n)];
dR = ([dq;du]-dys)/qscl;
end
%% Function for the calculation of qE at the next time level
function [qE,dqE] = solve_qE(qE,S,b,nonlinear_elements,q,u,Om,dt,...
        dS_dOm,dS_ddel,db,dq,du,dOm,ddel,fscl)

% Newton iteration with Cholesky factorization of Jacobian
for ik=1:10
    uE = ( 2/dt*(qE-q) - u )*Om;
    [fnl,dfnl_dq,dfnl_du] = nonlinear_forces(qE,uE,nonlinear_elements);
    dfnl_dqE = dfnl_dq+dfnl_du*2/dt*Om;
    Chol = chol(S+dfnl_dqE);
    qE = qE-Chol\(Chol'\(S*qE+fnl-b));
    r = (S*qE+fnl-b)/fscl;
    if max(abs(r))<1e-4
        break;
    end
end
if max(abs(r))>1e-4
    disp('noconv');
end

% Calculate gradient
dqE = Chol\(Chol'\(db+dfnl_du*...
    ( 2/dt*dq*Om + du*Om + (2/dt*(q-qE)+u)*dOm ) ...
    - dS_dOm*qE*dOm - dS_ddel*qE*ddel));
end
%% Function for the evaluation of the nonlinear forces
function [fnl,dfnl_dq,dfnl_du] = nonlinear_forces(q,u,nonlinear_elements)
% Initialize nonlinear force vector and its derivatives
fnl = zeros(size(q,1),1);
dfnl_dq = zeros(size(fnl,1),size(q,1));
dfnl_du = zeros(size(fnl,1),size(u,1));

% Evaluate all nonlinear elements
for nl=1:length(nonlinear_elements)
    % Determine force direction and calculate displacement and velocity of
    % nonlinear element
    w = nonlinear_elements{nl}.force_direction;
    qnl = w'*q; unl = w'*u;
    
    % Determine nonlinear force
    switch lower(nonlinear_elements{nl}.type)
        case 'cubicspring'
            fnl = fnl + ...
                w*nonlinear_elements{nl}.stiffness*qnl.^3;
            dfnl_dq = dfnl_dq + ...
                w*nonlinear_elements{nl}.stiffness*3*qnl.^2*w';
        case 'quadraticdamper'
            fnl = fnl + ...
                w*nonlinear_elements{nl}.damping*qnl.^2.*unl;
            dfnl_dq = dfnl_dq + ...
                w*nonlinear_elements{nl}.damping*2*qnl.*unl*w';
            dfnl_du = dfnl_du + ...
                w*nonlinear_elements{nl}.damping*qnl.^2*w';
        case 'unilateralspring'
            fnl = fnl + ...
                w*nonlinear_elements{nl}.stiffness*...
                (qnl-nonlinear_elements{nl}.gap).*...
                double(qnl-nonlinear_elements{nl}.gap>=0);
            dfnl_dq = dfnl_dq + ...
                w*nonlinear_elements{nl}.stiffness*...
                double(qnl-nonlinear_elements{nl}.gap>=0)*w';
        case 'tanhdryfriction'
            fnl = fnl + ...
                w*(nonlinear_elements{nl}.friction_limit_force.*...
                tanh(unl./nonlinear_elements{nl}.eps));
            dfnl_du = dfnl_du + ...
                w*(nonlinear_elements{nl}.friction_limit_force.*...
                ( 1 - tanh(unl./nonlinear_elements{nl}.eps).^2 )./...
                nonlinear_elements{nl}.eps)*w';
%         case 'elasticdryfriction': Hysteretic nonlinearities require
%         special treatment!
            
        case 'fe'
            myAssembly = nonlinear_elements{nl}.assembly;
            [dfnl_dq, fnl] = myAssembly.tangent_stiffness_and_force(...
                        myAssembly.unconstrain_vector(q) );
            dfnl_dq = myAssembly.constrain_matrix( dfnl_dq );
            fnl = myAssembly.constrain_vector( fnl );
            
            % remove linear terms
            K0 = nonlinear_elements{nl}.K;
            dfnl_dq = dfnl_dq - K0;
            fnl = fnl - K0*q;
            
        case 'custom'
            % myAssembly = nonlinear_elements{nl}.assembly;
            [dfnl_dq, fnl] = fnl_CUSTOM( q );
            
        otherwise
            error(['Unknown nonlinear element ' ...
                nonlinear_elements{nl}.type '.']);
    end
end
end