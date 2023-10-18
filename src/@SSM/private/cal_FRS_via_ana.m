function cal_FRS_via_ana(obj,parRange,gammas,lambda,r,outdof,optdof,optData,R1)
% This function calcualtes the forced response surface of the system via
% analytic prediction


optData.coordinates = 'polar';
%% zero function for FRS in reduced coordinates
isbaseForce = obj.System.Options.BaseExcitation;
frsFun = @(omega,epsilon,rho) frs_level_set(rho,omega,epsilon,gammas,lambda(1),abs(full(r)),isbaseForce);
%% plot FRS via fimplicit3
startfi  = tic;
rhomax   = obj.FRSOptions.rhoMax;
meshdens = obj.FRSOptions.meshDens;
figure;
fimplicit3(frsFun,[parRange{1}(1) parRange{1}(2) parRange{2}(1) parRange{2}(2) 0 rhomax], 'MeshDensity',meshdens);
xlabel('$\Omega$','FontSize',16,'Interpreter',"latex");
ylabel('$\epsilon$','FontSize',16,'Interpreter',"latex");
zlabel('$\rho$','FontSize',16,'Interpreter',"latex");
timings.fimp = toc(startfi);

% fimplicit does not return patches but it is robust to generate the FRS,
% from which we can tell what is a good esimation for rhomax. Then one can
% set up new rhomax for isosurface plot.

%% plot FRS via isosurface (zero contour)
oo = linspace(parRange{1}(1), parRange{1}(2),meshdens);
ee = linspace(parRange{2}(1), parRange{2}(2),meshdens);
rr = linspace(0,rhomax,meshdens);
[Om,Ep,Rho] = meshgrid(oo,ee,rr);
VV = frsFun(Om,Ep,Rho);
I  = isosurface(Om,Ep,Rho,VV,0);

figure;
patch('Faces', I.faces, 'Vertices', I.vertices, ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
view([20,50]); grid on; box on
xlabel('$\Omega$','FontSize',16,'Interpreter',"latex");
ylabel('$\epsilon$','FontSize',16,'Interpreter',"latex");
zlabel('$\rho$','FontSize',16,'Interpreter',"latex");
xlim([parRange{1}(1), parRange{1}(2)]);
ylim([parRange{2}(1), parRange{2}(2)]);


%% Map FRS from reduced coordinates to physical coordinates
% solutions in reduced coordinates
rho = I.vertices(:,3);
om  = I.vertices(:,1);
ep  = I.vertices(:,2);
[a,b] = frc_ab(rho,om,gammas,lambda(1));
th    = atan2(b.*real(r)-a.*imag(r), -a.*real(r)-b.*imag(r));
npts  = numel(rho);
isnonauto = obj.Options.contribNonAuto;
if isnonauto
    upo = [rho th om ep]; %[rho theta omega,epsilon]
else
    upo = [rho th om]; % [rho theta omega]
end

nout    = numel(outdof);
outidx  = zeros(nout,1);
for k=1:nout
    idxk = find(outdof(k)==optdof);
    assert(numel(idxk)>0,'outdof is not included in optdof');
    outidx(k) = idxk;
end
optData.outdof = outdof;
optData.outidx = outidx;
optData.mFreqs = 1;
optData.scales = ones(2*nout+1,1);

% calculate observables
yout = zeros(npts,2*nout+1);
for k = 1:npts
    uk = upo(k,:)';
    [~, y] = frs_output([], optData, uk);
    yout(k,:) = y';
end

figure;
patch('Faces', I.faces, 'Vertices', [om,ep,yout(:,1)], ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeAlpha', 0.4);
view([20,50]); grid on; box on
set(gca,'FontSize',12);
set(gca,'LineWidth',1.2);
xlabel('$\Omega$','FontSize',16,'Interpreter',"latex");
ylabel('$\epsilon$','FontSize',16,'Interpreter',"latex");
zlabel('$A_\mathrm{obj}$','FontSize',16,'Interpreter',"latex");
xlim([parRange{1}(1), parRange{1}(2)]);
ylim([parRange{2}(1), parRange{2}(2)]);


% plot FRS with stability
colors_chart = {'magenta',[0.7 0.7 0.7]};
nfaces = size(I.faces,1);
stab   = false(nfaces,1);
for k=1:nfaces
   vtx  = I.faces(k,:);
   rhok = rho(vtx);
   thk  = th(vtx);
   epk  = ep(vtx);
   stk  = check_stability(rhok,thk,gammas,lambda(1),epk,R1);
   stab(k) = all(stk);
end
figure; hold on
alidx  = 1:nfaces;
stidx  = find(stab);
ustidx = setdiff(alidx,stidx);
% stable part
patch('Faces', I.faces(stidx,:), 'Vertices', [om,ep,yout(:,1)], ...
  'FaceColor', colors_chart{2}, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)   
% unstable part
patch('Faces', I.faces(ustidx,:), 'Vertices', [om,ep,yout(:,1)], ...
  'FaceColor', colors_chart{1}, 'FaceAlpha', 0.4, 'EdgeAlpha', 0.2)   
hold off
view([20 50]); grid on; box on;
set(gca,'FontSize',12);
set(gca,'LineWidth',1.2);
xlabel('$\Omega$','FontSize',16,'Interpreter',"latex");
ylabel('$\epsilon$','FontSize',16,'Interpreter',"latex");
zlabel('$A_\mathrm{obj}$','FontSize',16,'Interpreter',"latex");
xlim([parRange{1}(1), parRange{1}(2)]);
ylim([parRange{2}(1), parRange{2}(2)]);

fprintf('the number of faces is %i\n',nfaces);

save('FRSana.mat','om','ep','rho','yout','outdof','stab','nfaces');

end