function y = reduced_dynamics_symbolic(lamdMaster,R0,options,varargin)
%REDUCED_DYNAMICS_SYMBOLIC This function returns vector field in symbolic
%expression, which can be used for latex documentation. In addition, it
%write the vector field to a function file in hard disk.
%
% y = reduced_dynamics_symbolic(lamdMaster,R0,type)
%
% lamdMaster - spectrum of master spectral subspace
% R0 - Reduced dynamics from compute_whisker
% options - a structure array with fields isauto,isdamped and numDigits. 
% isauto=true means that only autonomous part of reduced dynamics is considered.
% isdamped=true indicates that damping forces are included. You may set
% isauto=true and isdamped=false if you study backbone curves. numDigits
% represents the number of digits for coefficients

sympref('FloatingPointOutput',false);
lamdRe = real(lamdMaster);
lamdIm = imag(lamdMaster);
lamdRe = lamdRe(1:2:end);
lamdIm = lamdIm(1:2:end);
order  = numel(R0);
m      = numel(lamdRe);
beta   = cell(m,1); % coefficients - each cell corresponds to one mode
kappa  = cell(m,1); % exponants

% create symbolic array
rho = []; theta = [];
for k=1:m
    rho = [rho; sym(['rho_',num2str(k)])];
    theta = [theta; sym(['theta_',num2str(k)])];
end

% assemble and simplify nonlinear coefficients (when isdamped=false)
for k = 2:order
    R = R0{k};
    coeffs = R.coeffs;
    ind = R.ind;
    if ~isempty(coeffs)
        for i=1:m
            betai = coeffs(2*i-1,:);
            % remove small coefficients
            Rebetai = real(betai);
            Imbetai = imag(betai);
            if ~options.isdamped
                bounds  = max(1e-6*norm(betai),1e-8);
                Rebetai(abs(Rebetai)<bounds) = 0;
                Imbetai(abs(Imbetai)<bounds) = 0;
            end
            betai = Rebetai+1j*Imbetai;
            [~,ki,betai] = find(betai);
            kappai = ind(ki,:);
            % assemble terms
            beta{i}  = [beta{i} betai];
            kappa{i} = [kappa{i}; kappai];
        end
    end
end

% formualtion of vector field
% remove small parts of eigenvalues when isdamped=false
if ~options.isdamped
    lamdRe(abs(lamdRe)<1e-6*norm(lamdMaster)) = 0;
    lamdIm(abs(lamdIm)<1e-6*norm(lamdMaster)) = 0;
end
% autonomous part
yrho = lamdRe.*rho;
yth  = sym(lamdIm);
em = eye(m);
for i=1:m
    kappai = kappa{i};
    kappai = full(kappai);
    betai  = beta{i};
    nka = size(kappai,1);
    nbe = numel(betai);
    assert(nka==nbe, 'Size of kappa%d and beta%d does not match',i,i);
    ei = em(i,:);
    for k=1:nka
        ka = kappai(k,:);
        be = betai(k);
        l = ka(1:2:end-1);
        j = ka(2:2:end);
        ang = (l-j-ei)*theta;
        rhopower = rho.^((l+j)');
        pdrho = prod(rhopower,1);
        yrho(i) = yrho(i)+pdrho*(real(be)*cos(ang)-imag(be)*sin(ang));
        yth(i)  = yth(i)+pdrho/rho(i)*(real(be)*sin(ang)+imag(be)*cos(ang));
    end
end
% nonautonomous part
if ~options.isauto && ~isempty(varargin)
    syms Omega epsilon
    assume(epsilon,'real')
    % preprocessing nonautonomous part
    R1 = varargin{1}; mFreqs = varargin{2}; kappas = varargin{3};
    iNonauto = []; % indices for resonant happens
    rNonauto = []; % value of leading order contribution
    kNonauto = []; % (pos) kappa indices with resonance
    kappa_set= kappas; % each row corresponds to one kappa
    kappa_pos = kappa_set(kappa_set>0);
    num_kappa = numel(kappa_pos); % number of kappa pairs
    for k=1:num_kappa
        kappak = kappa_pos(k);
        idm = find(mFreqs(:)==kappak);
        R_10 = R1{1}.coeffs;
        idk = find(kappa_set==kappak);
        r = R_10(2*idm-1,idk);    
        iNonauto = [iNonauto; idm];
        rNonauto = [rNonauto; r];  
        kNonauto = [kNonauto; idk];
    end     
    % add contribution to vector field
    yth  = yth-mFreqs(:)*Omega;
    for i=1:numel(iNonauto)
        id = iNonauto(i);
        r  = epsilon*rNonauto(i);
        rRe = real(r);
        rIm = imag(r);
        yrho(id,:) = yrho(id,:)+rRe.*cos(theta(id,:))+rIm.*sin(theta(id,:));
        yth(id,:)  = yth(id,:)-rRe.*sin(theta(id,:))./rho(id,:)+rIm.*cos(theta(id,:))./rho(id,:);
    end      
end

% write to y
y = [yrho;yth];
y(1:2:end-1) = yrho;
y(2:2:end)   = yth;
y = vpa(y,options.numDigits);
end