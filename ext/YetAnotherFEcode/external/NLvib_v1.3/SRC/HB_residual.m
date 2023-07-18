%========================================================================
% DESCRIPTION: 
% Matlab function setting up the frequency-domain residual vector 'R' and 
% its derivatives for given frequency 'Om' and vector of harmonics of the 
% generalized coordiantes 'X'. The corresponding time-domain model 
% equation is 
% 
%       System.M * \ddot q + System.D * \dot q + System.K * q + ...
%                   f_nl(q,\dot q) - f_ex(t) = 0.
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
function [R,dR,Q] = ...
    HB_residual(X,system,H,N,analysis_type,varargin)
%% Handle input variables depending on the modus

% System matrices
M = system.M;
D = system.D;
K = system.K;

% number of degrees of freedom
n = size(M,1);

% Conversion from sine-cosine to complex-exponential representation
I0 = 1:n; ID = n+(1:H*n);
IC = n+repmat(1:n,1,H)+n*kron(0:2:2*(H-1),ones(1,n)); IS = IC+n;
dX = eye(length(X));
Q = zeros(n*(H+1),1);   dQ = zeros(size(Q,1),size(dX,2));
Q(I0) = X(I0);          dQ(I0,:) = dX(I0,:);
Q(ID) = X(IC)-1i*X(IS); dQ(ID,:) = dX(IC,:)-1i*dX(IS,:);

% Handle analysis type
if nargin<=4 || isempty(analysis_type)
    % Default analysis: frequency response
    analysis_type = 'frf';
end
switch lower(analysis_type)
    case {'frf','frequency response'}
        % Frequency response analysis: X = [Q;Om]
        
        % Excitation 'excitation' is the fundamental harmonic of the
        % external forcing
        Fex1 = system.Fex1;
        % Setup excitation vector
        Fex = zeros(n*(H+1),1);
        h = 1;% Forcing of h-th harmonic
        Fex(h*n+(1:n)) = Fex1;
        dFex = zeros(size(Fex,1),length(X));
        
        % Derivative of damping (does not depend on unknowns in the case of
        % frequency response)
        dD_dalpha = 0*M; dalpha = zeros(1,length(X));
        
        % Excitation frequency
        Om  =  X(end)	;
        dOm = dX(end,:)	;
        % Scaling of dynamic force equilibrium
        fscl = 1;
    case {'nma','nonlinear modal analysis'}
        % Nonlinear modal analysis:  X = [Psi;Om;del;log10a]
        
        % Interpret additional input
        inorm = varargin{1};
        
        % Modal mass
        a = exp(log(10)*X(end));
        da = log(10)*exp(log(10)*X(end))*dX(end,:);
        
        % In this case, the harmonics are amplitude-normalized. We thus
        % have to scale them by the amplitude a.
        Psi = Q; dPsi = dQ;
        Q = Psi*a;
        dQ = dPsi*a + Psi*da;
        
        % Modal frequency
        Om = X(end-2);
        dOm = dX(end-2,:);
        
        % Modal damping ratio
        del = X(end-1);
        ddel = dX(end-1,:);
        
        % Extended Periodic Motion concept: artifical negative viscous 
        % mass-proportional damping
        alpha = 2*Om*del;
        dalpha = 2*dOm*del+2*Om*ddel;
        D = D-alpha*M;
        dD_dalpha = -M;
        
        % No external forcing
        Fex = zeros(n*(H+1),1);
        dFex = zeros(size(Fex,1),length(X));
        
        % Scaling of dynamic force equilibrium
        if length(varargin)<2 || isempty(varargin{2})
            fscl = 1;
        else
            fscl = varargin{2};
        end
    otherwise
        error(['Unknown analysis type ' analysis.type '.']);
end
%% Computation of the Fourier coefficients of the nonlinear forces and the 
% Jacobian using AFT
[Fnl,dFnl] = HB_nonlinear_forces_AFT(Q,dQ,Om,dOm,H,N,...
    system.nonlinear_elements);
%% Assembly of the residual and the Jacobian

% Dynamic force equilibrium
Rc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*Q + Fnl - Fex;
dRc = ( -Om^2*kron(diag((0:H).^2),M) + 1i*Om*kron(diag(0:H),D) + ...
    kron(eye(H+1),K) )*dQ + dFnl - dFex + ...
    ( -2*Om*kron(diag((0:H).^2),M) + 1i*kron(diag(0:H),D) )*Q*dOm + ...
    1i*Om*kron(diag(0:H),dD_dalpha)*Q*dalpha;

% Scale dynamic force equilibrium (useful for numerical reasons)
Rc = 1/fscl*(Rc);
dRc = 1/fscl*(dRc);

% Conversion from complex-exponential to sine-cosine representation
R = zeros(size(X,1)-1,1); dR = zeros(size(X,1)-1,size(X,1));
R(I0) = real(Rc(I0)); dR(I0,:) = real(dRc(I0,:));
R(IC) = real(Rc(ID)); dR(IC,:) = real(dRc(ID,:));
R(IS) = -imag(Rc(ID)); dR(IS,:) = -imag(dRc(ID,:));

if strcmpi(analysis_type,'nma') || ...
        strcmpi(analysis_type,'nonlinear modal analysis')
    % Scale dynamic force equilibrium by modal amplitude
    % NOTE: We first evaluate the derivative, as we then overwrite R!
    dR(1:end-2,:) = dR(1:end-2,:)/a-R(1:end-2)/a^2*da;
    R(1:end-2) = R(1:end-2)/a;
    
    % Amplitude normalization: The mass of the nonlinear mode shape (all
    % harmonics) is enforced to be one.
    R(end-1) = real(Psi'*kron(eye(H+1),M)*Psi-1);
    dR(end-1,:) = real(2*(Psi'*kron(eye(H+1),M))*dPsi);
    
    % Phase normalization: Velocity of coordinate 'inorm' is enforced to be 
    % zero at t=0.
    R(end) = (1:H)*imag(Psi(inorm+(n:n:H*n)));
    dR(end,:) = (1:H)*imag(dPsi(inorm+(n:n:H*n),:));
end
end
%% Computation of the Fourier coefficients of the nonlinear forces and the 
% Jacobian using AFT
function [F,dF] = ...
    HB_nonlinear_forces_AFT(Q,dQ,Om,dOm,H,N,nonlinear_elements)
%% Initialize output
F = zeros(size(Q));
dF = zeros(size(F,1),size(dQ,2));

%% Iterate on nonlinear elements
for nl=1:length(nonlinear_elements)
    % Specify time samples along period
    tau = (0:2*pi/N:2*pi-2*pi/N)';
    
    % If the nonlinear element is hysteretic (e.g. elastic dry friction
    % element), we march 2 periods to reach a stablizized hysteresis.
    if nonlinear_elements{nl}.ishysteretic
        tau = [tau;tau+2*pi];
    end
    
    if nonlinear_elements{nl}.islocal
        % Standard case of a local nonlinearity describing a discrete
        % nonlinear element attached to a mechanical system, where the
        % contribution \Delta f to the global nonlinear force vector f is
        %   \Delta f = w * fnl(qnl,unl)
        % where the force fnl is scalar, acts in the direction w and 
        % depends on
        %   qnl=w'*q, unl = w'*u.
        
        % Determine force direction associated with nonlinear element
        if size(Q,1)==H+1
            w = 1;
        else
            w = nonlinear_elements{nl}.force_direction;
        end
        W = kron(eye(H+1),w);
        
        % Apply inverse discrete Fourier transform
        H_iDFT = exp(1i*tau*(0:H));
        qnl = real(H_iDFT*(W'*Q));
        dqnl = real(H_iDFT*(W'*dQ));
        
        %% Evaluate nonlinear force in time domain
        switch lower(nonlinear_elements{nl}.type)
            case 'cubicspring'
                fnl = nonlinear_elements{nl}.stiffness*qnl.^3;
                dfnl = nonlinear_elements{nl}.stiffness*3*...
                    repmat(qnl.^2,1,size(dqnl,2)).*dqnl;
            case 'quadraticdamper'
                % In this case, also the velocity is required
                % Frequency domain derivative matrix and its derivative w.r.t. Om
                qnldot = real(H_iDFT*(1i*Om*(0:H)'.*(W'*Q)));
                dqnldot = real(H_iDFT*(1i*Om*repmat((0:H)',1,size(dQ,2)).*...
                    (W'*dQ))) + real(H_iDFT*(1i*(0:H)'.*(W'*Q))*dOm);
                fnl = nonlinear_elements{nl}.damping*qnl.^2.*qnldot;
                dfnl = nonlinear_elements{nl}.damping*2*...
                    repmat(qnl.*qnldot,1,size(dqnl,2)).*dqnl + ...
                    nonlinear_elements{nl}.damping*...
                    repmat(qnl.^2,1,size(dqnl,2)).*dqnldot;
            case 'tanhdryfriction'
                % In this case, also the velocity is required
                % Frequency domain derivative matrix and its derivative w.r.t. Om
                qnldot = real(H_iDFT*(1i*Om*(0:H)'.*(W'*Q)));
                dqnldot = real(H_iDFT*(1i*Om*repmat((0:H)',1,size(dQ,2)).*...
                    (W'*dQ))) + real(H_iDFT*(1i*(0:H)'.*(W'*Q))*dOm);
                fnl = nonlinear_elements{nl}.friction_limit_force.*...
                    tanh(qnldot./nonlinear_elements{nl}.eps);
                dfnl = repmat(nonlinear_elements{nl}.friction_limit_force.*...
                    ( 1 - tanh(qnldot./nonlinear_elements{nl}.eps).^2 )./...
                    nonlinear_elements{nl}.eps,1,size(dqnldot,2)).*dqnldot;
            case 'unilateralspring'
                fnl = nonlinear_elements{nl}.stiffness*...
                    (qnl-nonlinear_elements{nl}.gap).*...
                    double(qnl-nonlinear_elements{nl}.gap>=0);
                dfnl = nonlinear_elements{nl}.stiffness*...
                    repmat(double(qnl-nonlinear_elements{nl}.gap>=0),...
                    1,size(dqnl,2)).*dqnl;
            case 'elasticdryfriction'
                % Set stiffness of elastic dry friction element (aka
                % Jenkins element)
                k = nonlinear_elements{nl}.stiffness;
                fc = nonlinear_elements{nl}.friction_limit_force;
                
                % Initialize force vector and Coulomb slider position
                fnl = zeros(size( qnl ));
                qsl = zeros(size( qnl ));
                dfnl = zeros(length(fnl),size(dQ,2));
                
                % Assume stuck conditions
                dfnl(1,:) = k*dqnl(1,:);
                
                % Iterative computation of force
                for ij = 2:length(fnl)
                    
                    % Predictor step
                    fnl(ij) = k*( qnl(ij) - qsl(ij-1) );
                    
                    % Corrector step
                    if(abs(fnl(ij)) >= fc)
                        fnl(ij) = fc*sign(fnl(ij));
                        qsl(ij) = qnl(ij) - fnl(ij)/k;
                    else
                        qsl(ij) = qsl(ij-1);
                        dfnl(ij,:) = k*(dqnl(ij,:)-dqnl(ij-1,:))+...
                            dfnl(ij-1,:);
                    end
                end
            otherwise
                error(['Unknown nonlinear element ' ...
                    nonlinear_elements{nl}.type '.']);
        end
        
        %% Forward Discrete Fourier Transform
        
        % Apply FFT
        Fnlc = fft(fnl(end-N+1:end))/N;
        dFnlc = fft(dfnl(end-N+1:end,:))/N;
        
        % Truncate and convert to half-spectrum notation
        Fnl = [real(Fnlc(1));2*Fnlc(2:H+1)];
        dFnl = [real(dFnlc(1,:));2*dFnlc(2:H+1,:)];
        
        % Store current force into global force vector
        F = F + W*Fnl;
        dF = dF + W*dFnl;
        
    else % Global nonlinearity
        switch lower(nonlinear_elements{nl}.type)
            case 'polynomialstiffness'
                % Inverse FFT
                Qc = transpose(reshape(Q,[],H+1)); n = size(Qc,2);
                Qc = [Qc(1,:); Qc(2:end,:)/2; zeros(N-H-1,n)];
                Qc(end-H+1:end,:) = flipud(conj(Qc(2:H+1,:)));
                q = ifft(Qc)*N;
                
                % Exponents and coefficients
                pp = nonlinear_elements{nl}.exponents;
                Et = transpose(nonlinear_elements{nl}.coefficients);
                
                % Evaluate polynimials
                z = prod(kron(q,ones(size(pp,1),1)) .^ repmat(pp,N,1),2);
                
                % Evaluate forces
                nz = size(Et,2);
                fnl = (Et*reshape(z,nz,N))';
                
                % Forward FFT
                Fnl = fft(fnl)/N;
                Fnl = [Fnl(1,:);Fnl(2:H+1,:)*2];
                
                % Add forces to global force vector
                F = F + reshape(transpose(Fnl),[],1);
                
                %% Analytical gradients 
                % (If the double loop bothers you, vectorize!)
                for l=1:n
                    % Apply IDFT to Jacobian of Q
                    ndx = size(dQ,2);
                    dQlc = [dQ(l,:); dQ(n+l:n:end,:)/2; zeros(N-H-1,ndx)];
                    dQlc(end-H+1:end,:) = flipud(conj(dQlc(2:H+1,:)));
                    dql = ifft(dQlc)*N;
                    
                    % Derive
                    notl = setdiff(1:n,l);
                    dzql_dql = repmat(pp(:,l),N,1).*...
                        (kron(q(:,l),ones(size(pp,1),1)).^repmat(pp(:,l)-1,N,1));
                    dzql_dql(isnan(dzql_dql)) = 0;
                    dz_dql = prod([ dzql_dql ...
                        kron(q(:,notl),ones(size(pp,1),1)).^...
                        repmat(pp(:,notl),N,1)], 2);
                    df_dql = (Et*reshape(dz_dql,nz,N))';
                    
                    % Derive
                    for i=1:n
                        dfi = repmat(df_dql(:,i),1,size(dql,2)).*dql;
                        dFi = fft(dfi(end-N+1:end,:))/N;
                        dFi = [dFi(1,:);dFi(2:H+1,:)*2];
                        dF(i:n:end,:) = dF(i:n:end,:) + dFi;
                    end
                end
                
            case 'fe'
                myAssembly = nonlinear_elements{nl}.assembly;
                K0 = nonlinear_elements{nl}.K;
                
                % Inverse FFT
                Qc = transpose(reshape(Q,[],H+1)); n = size(Qc,2);
                Qc = [Qc(1,:); Qc(2:end,:)/2; zeros(N-H-1,n)];
                Qc(end-H+1:end,:) = flipud(conj(Qc(2:H+1,:)));
                q = ifft(Qc)*N;
                
                % Analytical gradients ____________________________________
                % (If the double loop bothers you, vectorize!)
                fnl         = zeros(N,n);
                df_dql_all  = zeros(N,n,n);
                for iii = 1 : size( q , 1 )
                    x = q(iii, :)';                    
                    [Kt, fint] = myAssembly.tangent_stiffness_and_force(...
                        myAssembly.unconstrain_vector(x) );
                    df_dql_all(iii, :, :) = myAssembly.constrain_matrix(Kt) - K0;
                    fnl(iii, :) = myAssembly.constrain_vector(fint) - K0*x;                    
                end

                % Forward FFT
                Fnl = fft(fnl)/N;                   % harmonics * dofs = samples * dofs
                Fnl = [Fnl(1,:);Fnl(2:H+1,:)*2];    % selected harmonics * dofs

                % Add forces to global force vector
                F = F + reshape(transpose(Fnl),[],1);   % [F0_dof1, F0_dof2, ... , FH_dof(n-1), FH_dof(n)]

                for l_ind=1:n
                    % Apply IDFT to Jacobian of Q
                    ndx = size(dQ,2);
                    dQlc = [dQ(l_ind,:); dQ(n+l_ind:n:end,:)/2; zeros(N-H-1,ndx)];
                    dQlc(end-H+1:end,:) = flipud(conj(dQlc(2:H+1,:)));
                    dql = ifft(dQlc)*N;

                    df_dql = df_dql_all(:,:,l_ind);

                    % Derive
                    for i=1:n
                        dfi = repmat(df_dql(:,i),1,size(dql,2)).*dql;
                        dFi = fft(dfi(end-N+1:end,:))/N;
                        dFi = [dFi(1,:);dFi(2:H+1,:)*2];
                        dF(i:n:end,:) = dF(i:n:end,:) + dFi;
                    end
                end
                
            case 'custom'
                % myAssembly = nonlinear_elements{nl}.assembly;
                fnl_CUSTOM = nonlinear_elements{nl}.custom_function;
                
                % Inverse FFT
                Qc = transpose(reshape(Q,[],H+1)); n = size(Qc,2);
                Qc = [Qc(1,:); Qc(2:end,:)/2; zeros(N-H-1,n)];
                Qc(end-H+1:end,:) = flipud(conj(Qc(2:H+1,:)));
                q = ifft(Qc)*N;
                
                % Analytical gradients ____________________________________
                % (If the double loop bothers you, vectorize!)
                fnl         = zeros(N,n);
                df_dql_all  = zeros(N,n,n);
                for iii = 1 : size( q , 1 )
                    x = q(iii, :)';                    
                    [Kt, fint] = fnl_CUSTOM( x );
                    % note for future realese: fnl_CUSTOM( x, \dot x, \ddot x );
                    %                          (this will also require to
                    %                          change dF)
                    df_dql_all(iii, :, :) = Kt;
                    fnl(iii, :) = fint;                    
                end
                
                % Forward FFT
                Fnl = fft(fnl)/N;                   % harmonics * dofs = samples * dofs
                Fnl = [Fnl(1,:);Fnl(2:H+1,:)*2];    % selected harmonics * dofs

                % Add forces to global force vector
                F = F + reshape(transpose(Fnl),[],1);   % [F0_dof1, F0_dof2, ... , FH_dof(n-1), FH_dof(n)]

                for l_ind=1:n
                    % Apply IDFT to Jacobian of Q
                    ndx = size(dQ,2);
                    dQlc = [dQ(l_ind,:); dQ(n+l_ind:n:end,:)/2; zeros(N-H-1,ndx)];
                    dQlc(end-H+1:end,:) = flipud(conj(dQlc(2:H+1,:)));
                    dql = ifft(dQlc)*N;

                    df_dql = df_dql_all(:,:,l_ind);

                    % Derive
                    for i=1:n
                        dfi = repmat(df_dql(:,i),1,size(dql,2)).*dql;
                        dFi = fft(dfi(end-N+1:end,:))/N;
                        dFi = [dFi(1,:);dFi(2:H+1,:)*2];
                        dF(i:n:end,:) = dF(i:n:end,:) + dFi;
                    end
                end
                
            otherwise
                error(['Unknown global nonlinearity type ' ...
                    nonlinear_elements{nl}.type '.']);
        end
    end
end
end