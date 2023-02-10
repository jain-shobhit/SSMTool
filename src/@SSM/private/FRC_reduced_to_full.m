function FRC = FRC_reduced_to_full(obj,Nonauto,tb,FRC,FRCdata,W_0,W_1,outdof,varargin)
% FRC_REDUCED_TO_FULL This function maps the forced response curve of
% equilibria/periodic orbits/two-dimensional invariant tori to the forced
% response curve of periodic orbits/two/three-dimensional invariant tori of
% the full system

% pre-processing
nt       = obj.FRCOptions.nt;
iNonauto = Nonauto.iNonauto;
kNonauto = Nonauto.kNonauto;
rNonauto = Nonauto.rNonauto;
mFreqs   = FRCdata.mFreqs;
m        = numel(mFreqs);

% flag for saving ICs (used for numerical integration)
if numel(varargin)==2
    saveICflag = strcmp(varargin{2},'saveICs');
elseif numel(varargin)==1
    saveICflag = strcmp(varargin{1},'saveICs');
else
    saveICflag = false;
end
if saveICflag
    Z0_frc = []; % initial state
end
timeFRCPhysicsDomain = tic;

% check toolbox
if strcmp(tb, 'ep')
    isep = true;
    Zout_frc  = [];
    Znorm_frc = [];
    Aout_frc  = [];
else
    flag = strcmp(tb,'po') || strcmp(tb,'tor');
    assert(flag, 'toolbox should be ep/po/tor');
    isep = false;
    nlab = numel(FRC.lab);
    zTr  = cell(nlab,1);
    nSeg = FRC.nSeg;
    if isnumeric(outdof)
        noutdof = numel(outdof);
    elseif isa(outdof, 'function_handle') % function handle for observables
        noutdof = numel(outdof(zeros(size(W_0{1}.coeffs ,1),1)));
    end
end

% Loop around a resonant mode
om   = FRC.om;
epsf = FRC.ep;
if FRCdata.isomega && obj.Options.contribNonAuto
    if isempty(obj.E)
        obj.choose_E(FRCdata.modes);
    end
    lambda = obj.E.spectrum(1:2:end);
    lambda = lambda(mFreqs==1);
end
for j = 1:numel(om)
    % compute non-autonomous SSM coefficients
    obj.System.Omega = om(j);
    if FRCdata.isomega
        if obj.Options.contribNonAuto
            [W_1, R_1] = obj.compute_perturbed_whisker(FRCdata.order);

            R_10 = R_1{1}.coeffs;
            assert(~isempty(R_10), 'Near resonance does not occur, you may tune tol');
            f = R_10((kNonauto-1)*2*m+2*iNonauto-1);

            assert(norm(f-rNonauto)<1e-3*norm(f), 'inner resonance assumption does not hold');

            fprintf('the forcing frequency %.4d is nearly resonant with the eigenvalue %.4d + i%.4d\n', om(j), real(lambda(1)),imag(lambda(1)))
        else
            W_1 = [];
        end
    end
    % Forced response in Physical Coordinates
    if obj.System.Options.BaseExcitation
        epsf(j) = epsf(j)*(om(j))^2;
    end    
    %% ep toolbox
    if isep
        state = FRC.z(j,:);
        [Aout, Zout, z_norm, Zic] = compute_full_response_2mD_ReIm(W_0, W_1, state, epsf(j), nt, mFreqs, outdof);

        % collect output in array
        Aout_frc = [Aout_frc; Aout];
        Zout_frc = [Zout_frc; Zout];
        Znorm_frc = [Znorm_frc; z_norm];
    else   
    %% po toolbox and tor toolbox
        qTrj = FRC.qTr{j};
        tTrj = FRC.tTr{j};
        nt   = numel(tTrj);
        numSegs = nSeg(j);
        Zout_frc = zeros(nt,noutdof,numSegs);
        for k=1:numSegs
            xbp = qTrj(:,:,k);
            xbp = xbp';
            x_comp = xbp(1:2:end-1,:)+1i*xbp(2:2:end,:);
            state  = zeros(FRCdata.dim, nt); % state
            state(1:2:end-1,:) = x_comp;
            state(2:2:end,:)   = conj(x_comp);

            [~, Zout, ~,Zic] = compute_full_response_traj(W_0, W_1, epsf(j), tTrj, state, om(j), outdof);
            Zout_frc(:,:,k) = Zout';
        end
        zTr{j} = Zout_frc;

    end

    if saveICflag
        Z0_frc = [Z0_frc Zic]; % initial state
    end 
end
%% 
% Record output
if isep
    FRC.Aout_frc  = Aout_frc;
    FRC.Zout_frc  = Zout_frc;
    FRC.Znorm_frc = Znorm_frc;
else
    FRC.zTr  = zTr;
end
if saveICflag
    FRC.Z0_frc = Z0_frc; % initial state
end 
FRC.timeFRCPhysicsDomain = toc(timeFRCPhysicsDomain);
FRC.SSMorder   = FRCdata.order;
FRC.SSMnonAuto = obj.Options.contribNonAuto;
FRC.SSMispolar = FRCdata.ispolar;

end