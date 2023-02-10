function [B,A,fnl,fext,outdof] = build_model_1st(nElements1,nElements2)
%% Finite Element Set
% Geometry

l1 = 707; % mm
h1 = 1;   % mm
b1 = 10;   % mm 2
isViscoelastic = false;
% Mesh parameters

xOut1 = l1;
BC1 = 'C'; % cantilevered beam
[M1,C1,K1,fnl1,~,outdof1,nodes1,elements1] = beam_model(l1,b1,h1,BC1,nElements1,isViscoelastic,xOut1);
x1 = nodes1(:,1);
y1 = zeros(size(x1));
nodes1 = [x1, y1];
f1 = figure; PlotMesh(nodes1,elements1,0);

l2 = 1000; % mm
h2 = 1;
b2 = 1;  %2
xOut2 = l2;
BC2 = 'C'; % cantilevered beam
[M2,C2,K2,fnl2,~,outdof2,nodes2,elements2] = beam_model(l2,b2,h2,BC2,nElements2,isViscoelastic,xOut2);
y2 = nodes2(:,1) - l2;
x2 = l1 * ones(size(y2));
nodes2 = [x2 y2];
figure(f1); hold on; PlotMesh(nodes2,elements2,0);

% assembling constraints
n1 = size(M1,1);
n2 = size(M2,1);

nConstraints = 2;
constraints  = sparse(nConstraints,n1+n2);
constraints(1, outdof1(1)) = 1;  constraints(1, n1 + outdof2(2)) = 1;  % u1 = -w1
constraints(2, outdof1(2)) = 1;  constraints(2, n1 + outdof2(1)) = -1; % w1 = u1
G = constraints.';

M = blkdiag(M1,M2);   
C = blkdiag(C1,C2);
K = blkdiag(K1,K2);

zeroTmp1 = sparse(n1+n2,n1+n2);
zeroTmp2 = sparse(n1+n2,nConstraints);
zeroTmp3 = sparse(nConstraints,nConstraints);
A = [-K zeroTmp1 -G
    zeroTmp1 M zeroTmp2
    constraints zeroTmp2' zeroTmp3];
B = [C M zeroTmp2
    M zeroTmp1 zeroTmp2
    zeroTmp2' zeroTmp2' zeroTmp3]; 

% compute VMs
Mm = blkdiag(M1,M2,sparse(nConstraints,nConstraints));   
Km = [blkdiag(K1,K2), G;
    constraints, sparse(nConstraints,nConstraints)];
n_VMs = 10; % first n_VMs modes with lowest frequency calculated 
[V,omega2] = eigs(Km,Mm,n_VMs,'SM');
threshold = 1e16;
[V,omega2] = remove_stiff_modes(V,omega2,threshold);
[V,omega2] = remove_zero_modes(V,omega2);
omega = sqrt(diag(omega2));

% plot VMs
modes = [1 2]; 
for mode = modes
    u1 = [0 ; V(1:3:n1,mode)]; w1 = [0; V(2:3:n1,mode)];
    disp1 = [u1,w1];
    u2 = [0; V(n1+1:3:n1+n2,mode)]; w2 = [0; V(n1+2:3:n1+n2,mode)];
    disp2 = [-w2,u2];
    if mode==1
        factorMode = 5;
    else
        factorMode = 10;
    end
    figure;
    PlotFieldonDeformedMesh(nodes1,elements1,disp1,'factor',factorMode);
    hold on;
    PlotFieldonDeformedMesh(nodes2,elements2,disp2,'factor',factorMode);
    title(['Mode ' num2str(mode) ', Frequency =' num2str(omega(mode)/(2*pi)) ' Hz'])
end

%% assembling nonlinearity
n = 2*(n1+n2)+nConstraints;
fnl = cell(1,2);
for j = 1:length(fnl)
    sizej = n * ones(1,j+2);
    subsj = [fnl1{j}.subs; fnl2{j}.subs+n1];
    valsj = [fnl1{j}.vals; fnl2{j}.vals];    
    fnl{j} = sptensor(subsj,-valsj,sizej);
end

% Assembly forcing
% forcingDOF = outdof2(2);
xForcing = l1/2;
[~,forcenode] = min(abs(nodes1(:,1)- xForcing));
forcingDOF = (forcenode-1)*3-1; % transverse direction after applying cantilivered BC

% forcingDOF = outdof2(2);
f0 = zeros(n,1); 
f0(forcingDOF) = 1;
fext = f0;

% outdof = outdof1(1:2); 
% @SHOBHIT, PLEASE CHECK vertDOF below
vertpos = -l2/2;
[~,vertnode] = min(abs(nodes2(:,2)- vertpos));
vertDOF = n1+(vertnode-1)*3-1; % transverse direction (local)
outdof = [forcingDOF,vertDOF];
end



function [V,LAMBDA] = remove_stiff_modes(V,LAMBDA,lambdaThreshold)
%REMOVE_STIFF_MODES: This function removes modes with infinite eigenvalues
D = diag(LAMBDA);
D = abs(D);
n = numel(D);
idx0 = [find(isinf(D));find(isnan(D))]; % indices with inf and nan eigenvalues
D(idx0) = lambdaThreshold-1;
idx1 = find(D>lambdaThreshold);
idx2 = setdiff(1:n, [idx0;idx1]);
V = V(:,idx2);
LAMBDA = LAMBDA(idx2,idx2);
if ~isempty(idx0)
    fprintf('%i nan/inf eigenvalues are removed\n',numel(idx0));
end
if ~isempty(idx1)
    fprintf('%i eigenvalues with mangnitude larger than %d are removed\n',...
        numel(idx1), lambdaThreshold);
end
end

function [V,LAMBDA] = remove_zero_modes(V,LAMBDA)
%REMOVE_STIFF_MODES: This function removes modes with zero eigenvalues
D = diag(LAMBDA);
D = abs(D);
n = numel(D);
idx1 = find(D<eps*max(D)); % indices with zero eigenvalues
idx2 = setdiff(1:n, idx1);
V = V(:,idx2);
LAMBDA = LAMBDA(idx2,idx2);
if ~isempty(idx1)
    fprintf('%i zero eigenvalues are removed\n',numel(idx1));
end
end